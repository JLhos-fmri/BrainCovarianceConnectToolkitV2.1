function BCCT_CaSCN_Surf_compute_freesurfer(Parameter)

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
datainput = datainput';
%%
if GlobalScale % global mean scale
    DATINPUTTEMP = datainput;
    datainput = zeros(size(DATINPUTTEMP));
    for i = 1:size(datainput,2)
        DATtemp = zeros(1,size(datainput,1));
        DATtemp(find(MASK)) = DATINPUTTEMP(find(MASK),i);
        DATtemp = DATtemp/mean(DATtemp(find(MASK)));
        datainput(:,i) = DATtemp;
    end
    clear DATINPUTTEMP DATtemp;
end
SeedROItype = Parameter.SeedROItype;
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
RealCompPara.mod = 'surfmap';

MASK = logical(MASK);
maskedSignal = datainput(find(MASK),:)';
maskedsigname = fullfile(Outputdir,'maskedSignal.mat');
save(maskedsigname,'maskedSignal','MASK');
RealCompPara.MASK = MASK;
% Parameter.imageorder = get(D.ORDedit1,'string');
% Parameter.GCAorder = str2num(get(D.calordedit,'string'));
% Parameter.Calmethod1 = get(D.Calmethod1,'val');
% Parameter.perm = get(D.perm,'val');
% Parameter.permnum = str2num(get(D.permnum,'string'));

Imgorder = Parameter.imageorder;
if strcmp(Imgorder,'Defaults')
    IMGORD = 1:size(datainput,2);
else
    IMGORD = load(Imgorder);
end
RealCompPara.IMGORD = IMGORD;
GCAorder = Parameter.GCAorder;
RealCompPara.GCAorder = GCAorder;
PermMark = Parameter.perm;
PermNum = Parameter.permnum;
RealCompPara.PermMark = PermMark;
RealCompPara.PermNum = PermNum;
Calmethod = Parameter.Calmethod1;
RealCompPara.Calmethod1 = Calmethod;



[voll, Ml, mr_parmsl, volszl] = load_mgh([Indir,filesep,sublist{1},filesep,'surf',filesep,Marker{1}]);
[volr, Mr, mr_parmsr, volszr] = load_mgh([Indir,filesep,sublist{1},filesep,'surf',filesep,Marker{2}]);
N1size = length(voll);
N2size = length(volr);
RealCompPara.Ml = Ml;
RealCompPara.Mr = Mr;
RealCompPara.N1size = N1size;
RealCompPara.N2size = N2size;


save([Outputdir,filesep,'RealCompPara.mat'],'RealCompPara');

RealCompParasingle = RealCompPara;
RealCompParasingle.PermMark = 0;
for i = 1:size(ROIsignals,2)
    Ressingle(i) = BCCT_Perm_CaSCNMap_Cal(RealCompParasingle,maskedSignal,ROIsignals(:,i));
end

outmatfileGCA = [Outputdir,filesep,'TotalGCAsingle.mat'];
save(outmatfileGCA,'Ressingle');

if Calmethod
    for i = 1:size(ROIsignals,2)
        
        outnamefilex2y_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_X2Y.mgh'];
        outnamefiley2x_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_Y2X.mgh'];
        outnamefilex2yT_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_X2Y_trans.mgh'];
        outnamefiley2xT_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_Y2X_trans.mgh'];
        outnamefile5_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_NetX2Y.mgh'];
        outnamefilex2y_p_L = [Outputdir,filesep,'lh.Pval_ROI',sprintf('%05d',i),'_X2Y.mgh'];
        outnamefiley2x_p_L = [Outputdir,filesep,'lh.Pval_ROI',sprintf('%05d',i),'_Y2X.mgh'];
        
        outnamefilex2y_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_X2Y.mgh'];
        outnamefiley2x_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_Y2X.mgh'];
        outnamefilex2yT_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_X2Y_trans.mgh'];
        outnamefiley2xT_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_Y2X_trans.mgh'];
        outnamefile5_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_NetX2Y.mgh'];
        outnamefilex2y_p_R = [Outputdir,filesep,'rh.Pval_ROI',sprintf('%05d',i),'_X2Y.mgh'];
        outnamefiley2x_p_R = [Outputdir,filesep,'rh.Pval_ROI',sprintf('%05d',i),'_Y2X.mgh'];
        
        Dat = zeros(size(MASK));
        Dat(find(MASK)) = Ressingle(i).res.ResultMap1;
        save_mgh(Dat(1:N1size),outnamefilex2y_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2y_R,Mr);
%         SurfStatWriteData(outnamefilex2y,Dat);
            
        Dat = zeros(size(MASK));
        Dat(find(MASK)) = Ressingle(i).res.ResultMap2;
        save_mgh(Dat(1:N1size),outnamefiley2x_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2x_R,Mr);
%         SurfStatWriteData(outnamefiley2x,Dat);
        
        Dat = zeros(size(MASK));
        Dat(find(MASK)) = Ressingle(i).res.ResultMap3;
        save_mgh(Dat(1:N1size),outnamefilex2yT_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yT_R,Mr);
%         SurfStatWriteData(outnamefilex2yT,Dat);
        
        Dat = zeros(size(MASK));
        Dat(find(MASK)) = Ressingle(i).res.ResultMap4;
        save_mgh(Dat(1:N1size),outnamefiley2xT_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xT_R,Mr);
%         SurfStatWriteData(outnamefiley2xT,Dat);
        
        Dat = zeros(size(MASK));
        Dat(find(MASK)) = Ressingle(i).res.ResultMap5;
        save_mgh(Dat(1:N1size),outnamefile5_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),outnamefile5_R,Mr);
%         SurfStatWriteData(outnamefile5,Dat);
        
        Dat = zeros(size(MASK));
        Dat(find(MASK)) = Ressingle(i).res.ResultMap_p1;
        save_mgh(Dat(1:N1size),outnamefilex2y_p_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2y_p_R,Mr);
%         SurfStatWriteData(outnamefilex2y_p,Dat);
        
        Dat = zeros(size(MASK));
        Dat(find(MASK)) = Ressingle(i).res.ResultMap_p2;
        save_mgh(Dat(1:N1size),outnamefiley2x_p_L,Ml);
        save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2x_p_R,Mr);
%         SurfStatWriteData(outnamefiley2x_p,Dat);
        
        
    end
else
    for i = 1:size(ROIsignals,2)
        for i_ord = 1:GCAorder
            
            outnamefilex2y_L = [Outputdir,filesep,'lh.Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
            outnamefiley2x_L = [Outputdir,filesep,'lh.Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
            outnamefilex2yT_L = [Outputdir,filesep,'lh.Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
            outnamefiley2xT_L = [Outputdir,filesep,'lh.Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
            outnamefilex2yZ_L = [Outputdir,filesep,'lh.Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
            outnamefiley2xZ_L = [Outputdir,filesep,'lh.Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
            outnamefilex2yP_L = [Outputdir,filesep,'lh.Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
            outnamefiley2xP_L = [Outputdir,filesep,'lh.Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
            
            outnamefilex2y_R = [Outputdir,filesep,'rh.Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
            outnamefiley2x_R = [Outputdir,filesep,'rh.Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
            outnamefilex2yT_R = [Outputdir,filesep,'rh.Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
            outnamefiley2xT_R = [Outputdir,filesep,'rh.Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
            outnamefilex2yZ_R = [Outputdir,filesep,'rh.Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
            outnamefiley2xZ_R = [Outputdir,filesep,'rh.Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
            outnamefilex2yP_R = [Outputdir,filesep,'rh.Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
            outnamefiley2xP_R = [Outputdir,filesep,'rh.Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
            
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Ressingle(i).coef.ResultMap1{i_ord};
            save_mgh(Dat(1:N1size),outnamefilex2y_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2y_R,Mr);
%             SurfStatWriteData(outnamefilex2y,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Ressingle(i).coef.ResultMap2{i_ord};
            save_mgh(Dat(1:N1size),outnamefiley2x_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2x_R,Mr);
%             SurfStatWriteData(outnamefiley2x,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Ressingle(i).coef.res.T_T1{i_ord};
            save_mgh(Dat(1:N1size),outnamefilex2yT_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yT_R,Mr);
%             SurfStatWriteData(outnamefilex2yT,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Ressingle(i).coef.res.T_T2{i_ord};
            save_mgh(Dat(1:N1size),outnamefiley2xT_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xT_R,Mr);
%             SurfStatWriteData(outnamefiley2xT,Dat);
                       
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Ressingle(i).coef.res.Z_T1{i_ord};
            save_mgh(Dat(1:N1size),outnamefilex2yZ_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yZ_R,Mr);
%             SurfStatWriteData(outnamefilex2yZ,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Ressingle(i).coef.res.Z_T2{i_ord};
            save_mgh(Dat(1:N1size),outnamefiley2xZ_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xZ_R,Mr);
%             SurfStatWriteData(outnamefiley2xZ,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Ressingle(i).coef.res.P_T1{i_ord};
            save_mgh(Dat(1:N1size),outnamefilex2yP_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yP_R,Mr);
%             SurfStatWriteData(outnamefilex2yP,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Ressingle(i).coef.res.P_T2{i_ord};
            save_mgh(Dat(1:N1size),outnamefiley2xP_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xP_R,Mr);
%             SurfStatWriteData(outnamefiley2xP,Dat);
      
        end
    end
end

if PermMark
    disp('Now permutation start');
    for i = 1:size(ROIsignals,2)
        Res(i) = BCCT_Perm_CaSCNMap_Cal(RealCompPara,maskedSignal,ROIsignals(:,i));
    end
    
    outmatfileGCA = [Outputdir,filesep,'TotalGCA.mat'];
    save(outmatfileGCA,'Res');
    
    if Calmethod
        for i = 1:size(ROIsignals,2)
            
            outnamefilex2y_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_X2Y.mgh'];
            outnamefiley2x_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_Y2X.mgh'];
            outnamefilex2yT_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_X2Y_trans.mgh'];
            outnamefiley2xT_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_Y2X_trans.mgh'];
            outnamefile5_L = [Outputdir,filesep,'lh.ROI',sprintf('%05d',i),'_NetX2Y.mgh'];
            outnamefilex2y_p_L = [Outputdir,filesep,'lh.Pval_ROI',sprintf('%05d',i),'_X2Y.mgh'];
            outnamefiley2x_p_L = [Outputdir,filesep,'lh.Pval_ROI',sprintf('%05d',i),'_Y2X.mgh'];
            
            outnamefilex2y_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_X2Y.mgh'];
            outnamefiley2x_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_Y2X.mgh'];
            outnamefilex2yT_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_X2Y_trans.mgh'];
            outnamefiley2xT_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_Y2X_trans.mgh'];
            outnamefile5_R = [Outputdir,filesep,'rh.ROI',sprintf('%05d',i),'_NetX2Y.mgh'];
            outnamefilex2y_p_R = [Outputdir,filesep,'rh.Pval_ROI',sprintf('%05d',i),'_X2Y.mgh'];
            outnamefiley2x_p_R = [Outputdir,filesep,'rh.Pval_ROI',sprintf('%05d',i),'_Y2X.mgh'];
            
            if PermMark
                outnamefilex2y_perm_L = [Outputdir,filesep,'lh.Perm_ROI',sprintf('%05d',i),'_X2Y.mgh'];
                outnamefiley2x_perm_L = [Outputdir,filesep,'lh.Perm_ROI',sprintf('%05d',i),'_Y2X.mgh'];
                outnamefilex2yT_perm_L = [Outputdir,filesep,'lh.Perm_ROI',sprintf('%05d',i),'_X2Y_trans.mgh'];
                outnamefiley2xT_perm_L = [Outputdir,filesep,'lh.Perm_ROI',sprintf('%05d',i),'_Y2X_trans.mgh'];
                outnamefile5_perm_L = [Outputdir,filesep,'lh.Perm_ROI',sprintf('%05d',i),'_NetX2Y.mgh'];
                
                outnamefilex2y_perm_R = [Outputdir,filesep,'rh.Perm_ROI',sprintf('%05d',i),'_X2Y.mgh'];
                outnamefiley2x_perm_R = [Outputdir,filesep,'rh.Perm_ROI',sprintf('%05d',i),'_Y2X.mgh'];
                outnamefilex2yT_perm_R = [Outputdir,filesep,'rh.Perm_ROI',sprintf('%05d',i),'_X2Y_trans.mgh'];
                outnamefiley2xT_perm_R = [Outputdir,filesep,'rh.Perm_ROI',sprintf('%05d',i),'_Y2X_trans.mgh'];
                outnamefile5_perm_R = [Outputdir,filesep,'rh.Perm_ROI',sprintf('%05d',i),'_NetX2Y.mgh'];
            end
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Res(i).res.ResultMap1;
            save_mgh(Dat(1:N1size),outnamefilex2y_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2y_R,Mr);
            %         SurfStatWriteData(outnamefilex2y,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Res(i).res.ResultMap2;
            save_mgh(Dat(1:N1size),outnamefiley2x_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2x_R,Mr);
            %         SurfStatWriteData(outnamefiley2x,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Res(i).res.ResultMap3;
            save_mgh(Dat(1:N1size),outnamefilex2yT_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yT_R,Mr);
            %         SurfStatWriteData(outnamefilex2yT,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Res(i).res.ResultMap4;
            save_mgh(Dat(1:N1size),outnamefiley2xT_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xT_R,Mr);
            %         SurfStatWriteData(outnamefiley2xT,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Res(i).res.ResultMap5;
            save_mgh(Dat(1:N1size),outnamefile5_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefile5_R,Mr);
            %         SurfStatWriteData(outnamefile5,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Res(i).res.ResultMap_p1;
            save_mgh(Dat(1:N1size),outnamefilex2y_p_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2y_p_R,Mr);
            %         SurfStatWriteData(outnamefilex2y_p,Dat);
            
            Dat = zeros(size(MASK));
            Dat(find(MASK)) = Res(i).res.ResultMap_p2;
            save_mgh(Dat(1:N1size),outnamefiley2x_p_L,Ml);
            save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2x_p_R,Mr);
            %         SurfStatWriteData(outnamefiley2x_p,Dat);
            
            
            if PermMark
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).res.P_map1;
                save_mgh(Dat(1:N1size),outnamefilex2y_perm_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2y_perm_R,Mr);
                %             SurfStatWriteData(outnamefilex2y_perm,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).res.P_map2;
                save_mgh(Dat(1:N1size),outnamefiley2x_perm_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2x_perm_R,Mr);
                %             SurfStatWriteData(outnamefiley2x_perm,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).res.P_map3;
                save_mgh(Dat(1:N1size),outnamefilex2yT_perm_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yT_perm_R,Mr);
                %             SurfStatWriteData(outnamefilex2yT_perm,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).res.P_map4;
                save_mgh(Dat(1:N1size),outnamefiley2xT_perm_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xT_perm_R,Mr);
                %             SurfStatWriteData(outnamefiley2xT_perm,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).res.P_map5;
                save_mgh(Dat(1:N1size),outnamefile5_perm_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefile5_perm_R,Mr);
                %             SurfStatWriteData(outnamefile5_perm,Dat);
                
            end
        end
    else
        for i = 1:size(ROIsignals,2)
            for i_ord = 1:GCAorder
                
                outnamefilex2y_L = [Outputdir,filesep,'lh.Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                outnamefiley2x_L = [Outputdir,filesep,'lh.Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                outnamefilex2yT_L = [Outputdir,filesep,'lh.Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                outnamefiley2xT_L = [Outputdir,filesep,'lh.Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                outnamefilex2yZ_L = [Outputdir,filesep,'lh.Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                outnamefiley2xZ_L = [Outputdir,filesep,'lh.Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                outnamefilex2yP_L = [Outputdir,filesep,'lh.Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                outnamefiley2xP_L = [Outputdir,filesep,'lh.Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                
                outnamefilex2y_R = [Outputdir,filesep,'rh.Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                outnamefiley2x_R = [Outputdir,filesep,'rh.Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                outnamefilex2yT_R = [Outputdir,filesep,'rh.Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                outnamefiley2xT_R = [Outputdir,filesep,'rh.Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                outnamefilex2yZ_R = [Outputdir,filesep,'rh.Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                outnamefiley2xZ_R = [Outputdir,filesep,'rh.Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                outnamefilex2yP_R = [Outputdir,filesep,'rh.Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                outnamefiley2xP_R = [Outputdir,filesep,'rh.Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                if PermMark
                    outnamefilex2yPerm_L = [Outputdir,filesep,'lh.Coef_PermP_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                    outnamefiley2xPerm_L = [Outputdir,filesep,'lh.Coef_PermP_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                    
                    outnamefilex2yPerm_R = [Outputdir,filesep,'rh.Coef_PermP_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.mgh'];
                    outnamefiley2xPerm_R = [Outputdir,filesep,'rh.Coef_PermP_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.mgh'];
                end
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).coef.ResultMap1{i_ord};
                save_mgh(Dat(1:N1size),outnamefilex2y_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2y_R,Mr);
                %             SurfStatWriteData(outnamefilex2y,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).coef.ResultMap2{i_ord};
                save_mgh(Dat(1:N1size),outnamefiley2x_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2x_R,Mr);
                %             SurfStatWriteData(outnamefiley2x,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).coef.res.T_T1{i_ord};
                save_mgh(Dat(1:N1size),outnamefilex2yT_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yT_R,Mr);
                %             SurfStatWriteData(outnamefilex2yT,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).coef.res.T_T2{i_ord};
                save_mgh(Dat(1:N1size),outnamefiley2xT_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xT_R,Mr);
                %             SurfStatWriteData(outnamefiley2xT,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).coef.res.Z_T1{i_ord};
                save_mgh(Dat(1:N1size),outnamefilex2yZ_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yZ_R,Mr);
                %             SurfStatWriteData(outnamefilex2yZ,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).coef.res.Z_T2{i_ord};
                save_mgh(Dat(1:N1size),outnamefiley2xZ_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xZ_R,Mr);
                %             SurfStatWriteData(outnamefiley2xZ,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).coef.res.P_T1{i_ord};
                save_mgh(Dat(1:N1size),outnamefilex2yP_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yP_R,Mr);
                %             SurfStatWriteData(outnamefilex2yP,Dat);
                
                Dat = zeros(size(MASK));
                Dat(find(MASK)) = Res(i).coef.res.P_T2{i_ord};
                save_mgh(Dat(1:N1size),outnamefiley2xP_L,Ml);
                save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xP_R,Mr);
                %             SurfStatWriteData(outnamefiley2xP,Dat);
                
                if PermMark
                    Dat = zeros(size(MASK));
                    Dat(find(MASK)) = Res(i).coef.P_map1(i_ord,:);
                    save_mgh(Dat(1:N1size),outnamefilex2yPerm_L,Ml);
                    save_mgh(Dat(1+N1size:N1size+N2size),outnamefilex2yPerm_R,Mr);
                    %                 SurfStatWriteData(outnamefilex2yPerm,Dat);
                    
                    Dat = zeros(size(MASK));
                    Dat(find(MASK)) = Res(i).coef.P_map2(i_ord,:);
                    save_mgh(Dat(1:N1size),outnamefiley2xPerm_L,Ml);
                    save_mgh(Dat(1+N1size:N1size+N2size),outnamefiley2xPerm_R,Mr);
                    %                 SurfStatWriteData(outnamefiley2xPerm,Dat);
                    
                end
            end
        end
    end
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