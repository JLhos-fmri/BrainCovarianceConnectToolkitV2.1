function BCCT_MOD_Surf_compute_freesurfer(Parameter)

mfilenam = which('BCCT.m');
[pat nam ext] = fileparts(mfilenam);
Outputdir = Parameter.Outdir;
Indir = Parameter.Indir;
marker = Parameter.markerdir;
Marker = {marker,['r',marker(2:end)]};
FSlab = Parameter.FSlabs;
if FSlab(1)
    FStype = 'fsaverage';
elseif FSlab(2)
    FStype = 'fsaverage3';
elseif FSlab(3)
    FStype = 'fsaverage4';
elseif FSlab(4)
    FStype = 'fsaverage5';
elseif FSlab(5)
    FStype = 'fsaverage6';
elseif FSlab(6)
    FStype = 'fsaverage_sym';
end
[surfp,ab] = SurfStatReadSurf({[pat,filesep,'templates',filesep,FStype,filesep,'lh.pial'],[pat,filesep,'templates',filesep,FStype,filesep,'rh.pial']});
[surfw,ab] = SurfStatReadSurf({[pat,filesep,'templates',filesep,FStype,filesep,'lh.white'],[pat,filesep,'templates',filesep,FStype,filesep,'rh.white']});
surfshow = surfp;
surfshow.coord = (surfp.coord+surfw.coord)/2;
totalvertex = size(surfshow.coord,2);

GlobalScale = Parameter.GlobalScale;
masksdir = Parameter.masksdir;
if strcmp(masksdir,'Default')
    MASK = ones(1,totalvertex);
else
    MASK = load(masksdir);
    if length(MASK)~=totalvertex
        error('wrong mask selection')
    end
end
MASK(isnan(MASK)) = 0;
MASK(isinf(MASK)) = 0;
MASK = logical(MASK);
RealCompPara.MASK = MASK;
RealCompPara.mod = 'surfmap';
%%
folds = dir(Indir);
subnum = 1;
try
    for i = 1:length(folds)-2
        if folds(i+2).isdir % make sure it's directory
            foldtemp = folds(i+2).name;
            if length(foldtemp)>=9 % fsaverage length is 9
                if ~strcmp(foldtemp(1:9),'fsaverage') % it's not fsaverage or fsaverage*
                    if ~isempty(dir([Indir,filesep,foldtemp,filesep,'surf'])) % contain surf fold
                        sublist{subnum,1} = foldtemp;
                        subnum = subnum+1;
                    end
                end
            else
                if ~isempty(dir([Indir,filesep,foldtemp,filesep,'surf'])) % contain surf fold
                    sublist{subnum,1} = foldtemp;
                    subnum = subnum+1;
                end
            end
        end
    end
    for i = 1:length(sublist)
        datainput(i,:) = SurfStatReadData({[Indir,filesep,sublist{i},filesep,'surf',filesep,Marker{1}],[Indir,filesep,sublist{i},filesep,'surf',filesep,Marker{2}]});
    end
catch
    error('please check the input directory')
end
if size(datainput,2)~=totalvertex
    error('No match with the template');
end
datainput(isnan(datainput)) = 0;
datainput(isinf(datainput)) = 0;
datainput= datainput';
%%
if GlobalScale % global mean scale
    DATINPUTTEMP = datainput;
    datainput = zeros(size(DATINPUTTEMP));
    for i = 1:size(datainput,2)
        DATtemp = zeros(size(datainput,1),1);
        DATtemp(find(MASK)) = DATINPUTTEMP(find(MASK),i);
        DATtemp = DATtemp/mean(DATtemp(find(MASK)));
        datainput(:,i) = DATtemp;
    end
    clear DATINPUTTEMP DATtemp;
end
SeedROItype = Parameter.SeedROItype;
RealCompPara.SeedROItype = SeedROItype;
switch SeedROItype
    case 1 % MNI condition
        mnis = Parameter.SeedROI_mni;
        indexs = mniroi_fs(mnis,surfshow);
        ROIsignals = datainput(indexs,:)';
        Signalname = fullfile(Outputdir,'ROIsignal.mat');
        save(Signalname,'ROIsignals','indexs');
    case 2
        rois = Parameter.SeedROI_nifti;
        datroi = SurfStatReadData(rois);
        dataroi(isnan(dataroi)) = 0;
        dat00 = unique(dataroi);
        nrois = length(unique(dataroi))-1;
        for i = 1:nrois
            indexs = find(dataroi==dat00(i+1));
            ROIsignals(:,i) = mean(datainput(indexs,:),1)';
        end
        Signalname = fullfile(Outputdir,'ROIsignal.mat');
        save(Signalname,'ROIsignals');
    case 3
        mats = Parameter.SeedROI_mat;
        temps = load(Parameter.SeedROI_mat);
        varname = Parameter.SeedROI_varname;
        eval(['ROIsignals = temps.',varname,';']);
        Signalname = fullfile(Outputdir,'ROIsignal.mat');
        save(Signalname,'ROIsignals');
    case 4
        txts = Parameter.SeedROI_txt;
        ROIsignals = load(txts);
        Signalname = fullfile(Outputdir,'ROIsignal.mat');
        save(Signalname,'ROIsignals');
    otherwise
        disp('error ROI defined')
end

COVcond = Parameter.covlab;
RealCompPara.COVcond = COVcond;
COVdir = Parameter.COVdir;
RealCompPara.COVdir = COVdir;
if COVcond(1)==1
    COV = load(COVdir);
    RealCompPara.COVs = COV;
end

Modufactor = load(Parameter.ModuFactor);
if COVcond(1)==1
    COVnoInt = term(COV);
end
RealCompPara.Modufactor = Modufactor;
MODFACTOR = term(Modufactor);


[voll, Ml, mr_parmsl, volszl] = load_mgh([Indir,filesep,sublist{1},filesep,'surf',filesep,Marker{1}]);
[volr, Mr, mr_parmsr, volszr] = load_mgh([Indir,filesep,sublist{1},filesep,'surf',filesep,Marker{2}]);
N1size = length(voll);
N2size = length(volr);
RealCompPara.Ml = Ml;
RealCompPara.Mr = Mr;
RealCompPara.N1size = N1size;
RealCompPara.N2size = N2size;

save([Outputdir,filesep,'RealCompPara.mat'],'RealCompPara')

MASK = logical(MASK);
maskedSignal = datainput(find(MASK),:)';
maskedsigname = fullfile(Outputdir,'maskedSignal.mat');
save(maskedsigname,'maskedSignal','MASK');

for i = 1:size(ROIsignals,2)
    if i<10
        Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI00000',num2str(i),'.mgh']);
        OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI00000',num2str(i),'.mgh']);
        OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI00000',num2str(i),'.mgh']);
        
        Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI00000',num2str(i),'.mgh']);
        OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI00000',num2str(i),'.mgh']);
        OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI00000',num2str(i),'.mgh']);
        OutfilenametempMat = fullfile(Outputdir,['mat_ROI00000',num2str(i),'.mat']);
    elseif i<100
        Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI0000',num2str(i),'.mgh']);
        OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI0000',num2str(i),'.mgh']);
        OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI0000',num2str(i),'.mgh']);
        OutfilenametempMat = fullfile(Outputdir,['mat_ROI0000',num2str(i),'.mat']);
        Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI0000',num2str(i),'.mgh']);
        OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI0000',num2str(i),'.mgh']);
        OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI0000',num2str(i),'.mgh']);
    elseif i<1000
        Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI000',num2str(i),'.mgh']);
        OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI000',num2str(i),'.mgh']);
        OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI000',num2str(i),'.mgh']);
        OutfilenametempMat = fullfile(Outputdir,['mat_ROI000',num2str(i),'.mat']);
        Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI000',num2str(i),'.mgh']);
        OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI000',num2str(i),'.mgh']);
        OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI000',num2str(i),'.mgh']);
    else
        Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI00',num2str(i),'.mgh']);
        OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI00',num2str(i),'.mgh']);
        OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI00',num2str(i),'.mgh']);
        OutfilenametempMat = fullfile(Outputdir,['mat_ROI00',num2str(i),'.mat']);
        Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI00',num2str(i),'.mgh']);
        OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI00',num2str(i),'.mgh']);
        OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI00',num2str(i),'.mgh']);
    end
    
    
    SEED = term(ROIsignals(:,i));
    if COVcond(1)==1
        mod = 1+SEED+MODFACTOR+SEED*MODFACTOR+COVnoInt;
        contrast = [0,0,0.1,zeros(1,size(COVnoInt,2))];
    else
        mod = 1+SEED+MODFACTOR+SEED*MODFACTOR;
        contrast = [0,0,0,1];
    end
    for j = 1:size(maskedSignal,2)
        Y = maskedSignal(:,j);
        slm = SurfStatLinMod( Y, mod );
        slmt = SurfStatT(slm,contrast);
        T(j) = slmt.t;
        [Z(j) P(j)] = AS_TFRtoZ(slmt.t,'T',slmt.df,[]);
    end
    
    Dat = zeros(size(MASK));
    Dat(find(MASK)) = T;
        save_mgh(Dat(1:N1size),Outfilenametemp_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),Outfilenametemp_R,Mr);
%     SurfStatWriteData(Outfilenametemp,Dat,'a')
    Dat = zeros(size(MASK));
    Dat(find(MASK)) = P;
        save_mgh(Dat(1:N1size),OutfilenametempP_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),OutfilenametempP_R,Mr);
%     SurfStatWriteData(OutfilenametempP,Dat,'a')
    Dat = zeros(size(MASK));
    Dat(find(MASK)) = Z;
        save_mgh(Dat(1:N1size),OutfilenametempZ_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),OutfilenametempZ_R,Mr);
%     SurfStatWriteData(OutfilenametempZ,Dat,'a')
    save(OutfilenametempMat,'T','P','Z','MASK');
end
uiwait(msgbox('Finished'));
end

function ind = mniroi_fs(mnis,surfshow)
MNIS = zeros(size(surfshow.coord));
MNIS(1,:) = mnis(1);
MNIS(2,:) = mnis(2);
MNIS(3,:) = mnis(3);
dis = sqrt((MNIS(1,:)-surfshow.coord(1,:)).^2+(MNIS(2,:)-surfshow.coord(2,:)).^2+(MNIS(3,:)-surfshow.coord(3,:)).^2);
[mindis minind] = min(dis);
ind = minind(1);
end


