function BCCT_Mod_matrixcomputemain(Parameter)
pathfiles = which('BCCT.m');
[path nam ext] = fileparts(pathfiles);
Outputdir = Parameter.Outputdir;
RealCompPara.Outputdir = Outputdir;
Inputdir = Parameter.Inputdir;
RealCompPara.Inputdir = Inputdir;
SeedROItype = Parameter.SeedROItype;
RealCompPara.SeedROItype = SeedROItype;
Partinfo = Parameter.Partinfo;
RealCompPara.Partinfo = Partinfo;
GlobalScale = Parameter.GlobalScal;
RealCompPara.GlobaScale = GlobalScale;
RealCompPara.mod = 'volroi';

if SeedROItype<3
    [vinput,datainput] = Dynamic_read_dir_NIFTI(Inputdir);
    dims = vinput.dim;
    Ttrans = vinput.mat;
    VOXELSIZE = abs([vinput.mat(1,1),vinput.mat(2,2),vinput.mat(3,3)]);
    if strcmp(Parameter.masks,'Defaults')
        maskfiles = fullfile(path,'templates','grey.nii');
        [vmask,datamask] = Dynamic_read_dir_NIFTI(maskfiles);
        datamask = datamask>0.2;
    else
        try
            maskfiles = Parameter.masks;
            [vmask,datamask] = Dynamic_read_dir_NIFTI(maskfiles);
        catch
            error('Reselect the mask or fill Defaults');
        end
    end
    if any(vmask.dim-vinput(1).dim)
        Maskimg0 = fullfile(Outputdir,'mask.nii');
        %     (PI,PO,NewVoxSize,hld,TargetSpace)
        dynamicBC_Reslice(maskfiles,Maskimg0,vinput(1).dim,0,vinput(1).fname);
        clear vmask datamask
        [vmask,datamask] = Dynamic_read_dir_NIFTI(Maskimg0);
        if strcmp(Parameter.masks,'Defaults')
            datamask = datamask>0.2;
        end
        MaskimgU = fullfile(Outputdir,'maskused.nii');
        DATMASK = reshape(datamask,dims);
        DynamicBC_write_NIFTI(DATMASK,vmask,MaskimgU);
    end
    DATMASK = reshape(datamask,dims);
    DATMASK = logical(DATMASK);
    if GlobalScale
        DATINPUTTEMP = datainput;
        datainput = zeros(size(DATINPUTTEMP));
        for i = 1:size(datainput,2)
            DATtemp = zeros(size(datainput,1),1);
            DATtemp(find(DATMASK)) = DATINPUTTEMP(find(DATMASK),i);
            DATtemp = DATtemp/mean(DATtemp(find(DATMASK)));
            datainput(:,i) = DATtemp;
        end
        clear DATINPUTTEMP DATtemp;
    end
    dims = vinput.dim;
    Ttrans = vinput.mat;
    VOXELSIZE = abs([vinput(1).mat(1,1),vinput(1).mat(2,2),vinput(1).mat(3,3)]);
end


switch SeedROItype
    case 1 % MNI condition
        roipath = Parameter.SeedROI_seproi;
        [vroi,dataroi,roinamelist] = Dynamic_read_dir_NIFTI(roipath);
        if any(vroi.dim-vinput(1).dim)
            error('wrong dimension for the ROI defination!');
        end
        nrois = size(dataroi,2);
        for i = 1:nrois
            indexs = find(dataroi(:,i));
            ROIsignals(:,i) = mean(datainput(indexs,:),1)';
        end
        Signalname = fullfile(Outputdir,'ROIsignal.mat');
        Signalnamelist = fullfile(Outputdir,'roinamelist.mat');
        save(Signalname,'ROIsignals');
        save(Signalnamelist,'roinamelist');
    case 2
        rois = Parameter.SeedROI_nifti;
        [vroi,dataroi] = Dynamic_read_dir_NIFTI(rois);
        if any(vroi.dim-vinput(1).dim)
            Roisimg0 = fullfile(Outputdir,'roi.nii');
            dynamicBC_Reslice(rois,Roisimg0,vinput.dim,0,vinput.fname);
            [vroi,dataroi] = Dynamic_read_dir_NIFTI(Roisimg0);
        end
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
COVcond = Parameter.covs;
RealCompPara.COVcond = COVcond;
COVdir = Parameter.COVtext;
RealCompPara.COVdir = COVdir;
if COVcond==1
    COV = load(COVdir);
    RealCompPara.COVs = COV;
end
Modufactor = load(Parameter.ModuFactor);
RealCompPara.Modufactor = Modufactor;
save(fullfile(Outputdir,'RealCompPara.mat'),'RealCompPara');


%%
MODFACTOR = term(Modufactor);
if COVcond==1
    COVnoInt = term(COV);
end
for i = 1:size(ROIsignals,2)
    for j = 1:size(ROIsignals,2)
        SEED = term(ROIsignals(:,i));
        if COVcond==1
            mod = 1+SEED+MODFACTOR+SEED*MODFACTOR+COVnoInt;
            contrast = [0,0,0.1,zeros(1,size(COVnoInt,2))];
        else
            mod = 1+SEED+MODFACTOR+SEED*MODFACTOR;
            contrast = [0,0,0,1];
        end
        Y = ROIsignals(:,j);
        slm = SurfStatLinMod( Y, mod );
        slmt = SurfStatT(slm,contrast);
        T(i,j) = slmt.t;
        [Z(i,j) P(i,j)] = AS_TFRtoZ(slmt.t,'T',slmt.df,[]);
    end
end

Outputfile = fullfile(Outputdir,'MOD_tzp_res.mat');
save(Outputfile,'T','P','Z', 'slmt')
disp('Compute Finished!')
end

