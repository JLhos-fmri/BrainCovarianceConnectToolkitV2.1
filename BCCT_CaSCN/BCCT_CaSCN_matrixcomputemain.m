function BCCT_CaSCN_matrixcomputemain(Parameter)
pathfiles = which('BCCT.m');
[path nam ext] = fileparts(pathfiles);
Outputdir = Parameter.Outputdir;
RealCompPara.Outputdir = Outputdir;
Inputdir = Parameter.Inputdir;
RealCompPara.Inputdir = Inputdir;
SeedROItype = Parameter.SeedROItype;
RealCompPara.SeedROItype = SeedROItype;
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
save(fullfile(Outputdir,'RealCompPara.mat'),'RealCompPara');
% ROIsignals = zscore(ROIsignals);
for i = 1:size(ROIsignals,2)
    Res(i) = BCCT_Perm_CaSCNMap_Cal(RealCompPara,ROIsignals,ROIsignals(:,i));
end
% save test
if RealCompPara.Calmethod1
    for i = 1:size(ROIsignals,2)
        GCAres.GCA_res_x2y(i,:) = Res(i).res.ResultMap1; %x2y
        GCAres.GCA_res_y2x(i,:) = Res(i).res.ResultMap2; %y2x
        GCAres.GCA_res_x2y_trans(i,:) = Res(i).res.ResultMap3; %x2y_trans
        GCAres.GCA_res_y2x_trans(i,:) = Res(i).res.ResultMap4; %y2x_trans
        GCAres.GCA_res_fx2y(i,:) = Res(i).res.ResultMap5; %Fx2y
        GCAres.GCA_res_x2y_pval(i,:) = Res(i).res.ResultMap_p1;
        GCAres.GCA_res_y2x_pval(i,:) = Res(i).res.ResultMap_p2;
        if PermMark
            GCAres.GCA_res_x2y_permpval(i,:) = Res(i).res.P_map1;
            GCAres.GCA_res_y2x_permpval(i,:) = Res(i).res.P_map2;
            GCAres.GCA_res_x2y_trans_permpval(i,:) = Res(i).res.P_map3;
            GCAres.GCA_res_y2x_trans_permpval(i,:) = Res(i).res.P_map4;
            GCAres.GCA_res_fx2y_permpval(i,:) = Res(i).res.P_map5;
        end
    end
else
    for iord = 1:GCAorder
        for i = 1:size(ROIsignals,2)
            
            GCAres(iord).GCA_coef_x2y(i,:) = Res(i).coef.ResultMap1{iord}; %x2y
            GCAres(iord).GCA_coef_y2x(i,:) = Res(i).coef.ResultMap2{iord}; %y2x
            GCAres(iord).GCA_coef_x2yT(i,:) = Res(i).coef.res.T_T1{iord};
            GCAres(iord).GCA_coef_y2xT(i,:) = Res(i).coef.res.T_T2{iord};
            GCAres(iord).GCA_coef_x2yZ(i,:) = Res(i).coef.res.Z_T1{iord};
            GCAres(iord).GCA_coef_y2xZ(i,:) = Res(i).coef.res.Z_T2{iord};
            GCAres(iord).GCA_coef_x2y_pval(i,:) = Res(i).coef.res.P_T1{iord};
            GCAres(iord).GCA_coef_y2x_pval(i,:) = Res(i).coef.res.P_T2{iord};
            if PermMark
                GCAres(iord).GCA_coef_x2y_permpval(i,:) = Res(i).coef.P_map1(iord,:);
                GCAres(iord).GCA_coef_y2x_permpval(i,:) = Res(i).coef.P_map2(iord,:);
            end
        end
    end
    
end

save(fullfile(Outputdir,'ResultGCA'),'GCAres')

disp('Compute Finished!')
end