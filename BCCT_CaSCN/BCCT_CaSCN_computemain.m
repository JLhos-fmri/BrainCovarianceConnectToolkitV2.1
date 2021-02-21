function BCCT_CaSCN_computemain(Parameter)
pathfiles = which('BCCT.m');
[path nam ext] = fileparts(pathfiles);
Outputdir = Parameter.Outputdir;
RealCompPara.Outputdir = Outputdir;
Inputdir = Parameter.Inputdir;
RealCompPara.Inputdir = Inputdir;
SeedROItype = Parameter.SeedROItype;
RealCompPara.SeedROItype = SeedROItype;
GlobalScale = Parameter.GlobalScal;
RealCompPara.GlobalScale = GlobalScale;
RealCompPara.mod = 'volmap';

[vinput,datainput] = Dynamic_read_dir_NIFTI(Inputdir);
dims = vinput(1).dim;
Ttrans = vinput(1).mat;
VOXELSIZE = abs([vinput(1).mat(1,1),vinput(1).mat(2,2),vinput(1).mat(3,3)]);
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
RealCompPara.V = vmask;
% RealCompPara.DATMASK = DATMASK;
%% add in 2020-11-17
if GlobalScale % global mean scale
    DATINPUTTEMP = datainput;
    datainput = zeros(size(DATINPUTTEMP));
    for i = 1:size(datainput,2)
        DATtemp = zeros(size(datainput,1),1);
        DATtemp(find(DATMASK)) = DATINPUTTEMP(find(DATMASK),i);
        DATtemp = DATtemp/mean(find(DATMASK));
        datainput(:,i) = DATtemp;
    end    
    clear DATINPUTTEMP DATtemp;
end
%%

switch SeedROItype
    case 1 % MNI condition
        mnis = Parameter.SeedROI_mni;
        radius = Parameter.SeedROI_radius;
        indexs = mniroi(mnis,radius,Ttrans,dims,VOXELSIZE,DATMASK);
        ROIsignals = mean(datainput(indexs,:),1)';
        DATROI = zeros(dims);
        DATROI(indexs) = 1;
        ROIname = fullfile(Outputdir,'roi.nii');
        DynamicBC_write_NIFTI(DATROI,vmask,ROIname);
        Signalname = fullfile(Outputdir,'ROIsignal.mat');
        save(Signalname,'ROIsignals');
    case 2
        rois = Parameter.SeedROI_nifti;
        [vroi,dataroi] = Dynamic_read_dir_NIFTI(rois);
        if any(vroi.dim-vinput.dim)
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

DATMASK = logical(DATMASK);
maskedSignal = datainput(find(DATMASK),:)';
maskedsigname = fullfile(Outputdir,'maskedSignal.mat');
save(maskedsigname,'maskedSignal','DATMASK');

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
save([Outputdir,filesep,'RealCompPara.mat'],'RealCompPara');
RealCompParasingle = RealCompPara;
RealCompParasingle.PermMark = 0;
for i = 1:size(ROIsignals,2)
    Ressingle(i) = BCCT_Perm_CaSCNMap_Cal(RealCompParasingle,maskedSignal,ROIsignals(:,i));
end

outmatfileGCA = [Outputdir,filesep,'TotalGCAsingle.mat'];
save(outmatfileGCA,'Ressingle');
% save(Ressingle 'Ressingle')
if Calmethod
    for i = 1:size(ROIsignals,2)
        outnamefilex2y = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_X2Y.nii'];
        outnamefiley2x = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_Y2X.nii'];
        outnamefilex2yT = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_X2Y_trans.nii'];
        outnamefiley2xT = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_Y2X_trans.nii'];
        outnamefile5 = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_NetX2Y.nii'];
        outnamefilex2y_p = [Outputdir,filesep,'Pval_ROI',sprintf('%05d',i),'_X2Y.nii'];
        outnamefiley2x_p = [Outputdir,filesep,'Pval_ROI',sprintf('%05d',i),'_Y2X.nii'];
     
        Dat = zeros(vinput(1).dim);
        Dat(find(DATMASK)) = Ressingle(i).res.ResultMap1;
        DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2y);
        Dat = zeros(vinput(1).dim);
        Dat(find(DATMASK)) = Ressingle(i).res.ResultMap2;
        DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2x);
        Dat = zeros(vinput(1).dim);
        Dat(find(DATMASK)) = Ressingle(i).res.ResultMap3;
        DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yT);
        Dat = zeros(vinput(1).dim);
        Dat(find(DATMASK)) = Ressingle(i).res.ResultMap4;
        DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xT);
        Dat = zeros(vinput(1).dim);
        Dat(find(DATMASK)) = Ressingle(i).res.ResultMap5;
        DynamicBC_write_NIFTI(Dat,vinput(1),outnamefile5);
        Dat = zeros(vinput(1).dim);
        Dat(find(DATMASK)) = Ressingle(i).res.ResultMap_p1;
        DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2y_p);
        Dat = zeros(vinput(1).dim);
        Dat(find(DATMASK)) = Ressingle(i).res.ResultMap_p2;
        DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2x_p);
    end
else
    for i = 1:size(ROIsignals,2)
        for i_ord = 1:RealCompPara.GCAorder
            
            outnamefilex2y = [Outputdir,filesep,'Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.nii'];
            outnamefiley2x = [Outputdir,filesep,'Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.nii'];
            outnamefilex2yT = [Outputdir,filesep,'Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.nii'];
            outnamefiley2xT = [Outputdir,filesep,'Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.nii'];
            outnamefilex2yZ = [Outputdir,filesep,'Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.nii'];
            outnamefiley2xZ = [Outputdir,filesep,'Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.nii'];
            outnamefilex2yP = [Outputdir,filesep,'Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.nii'];
            outnamefiley2xP = [Outputdir,filesep,'Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.nii'];
            
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Ressingle(i).coef.ResultMap1{i_ord};
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2y);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Ressingle(i).coef.ResultMap2{i_ord};
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2x);
            
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Ressingle(i).coef.res.T_T1{i_ord};
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yT);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Ressingle(i).coef.res.T_T2{i_ord};
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xT);
            
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Ressingle(i).coef.res.Z_T1{i_ord};
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yZ);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Ressingle(i).coef.res.Z_T2{i_ord};
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xZ);
            
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Ressingle(i).coef.res.P_T1{i_ord};
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yP);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Ressingle(i).coef.res.P_T2{i_ord};
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xP);
        end
    end
end

if RealCompPara.PermMark
    disp('Now permutation start');
    for i = 1:size(ROIsignals,2)
        Res(i) = BCCT_Perm_CaSCNMap_Cal(RealCompPara,maskedSignal,ROIsignals(:,i));
    end
    if Calmethod
        for i = 1:size(ROIsignals,2)
            outnamefilex2y = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_X2Y.nii'];
            outnamefiley2x = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_Y2X.nii'];
            outnamefilex2yT = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_X2Y_trans.nii'];
            outnamefiley2xT = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_Y2X_trans.nii'];
            outnamefile5 = [Outputdir,filesep,'ROI',sprintf('%05d',i),'_NetX2Y.nii'];
            outnamefilex2y_p = [Outputdir,filesep,'Pval_ROI',sprintf('%05d',i),'_X2Y.nii'];
            outnamefiley2x_p = [Outputdir,filesep,'Pval_ROI',sprintf('%05d',i),'_Y2X.nii'];
            if PermMark
                outnamefilex2y_perm = [Outputdir,filesep,'Perm_ROI',sprintf('%05d',i),'_X2Y.nii'];
                outnamefiley2x_perm = [Outputdir,filesep,'Perm_ROI',sprintf('%05d',i),'_Y2X.nii'];
                outnamefilex2yT_perm = [Outputdir,filesep,'Perm_ROI',sprintf('%05d',i),'_X2Y_trans.nii'];
                outnamefiley2xT_perm = [Outputdir,filesep,'Perm_ROI',sprintf('%05d',i),'_Y2X_trans.nii'];
                outnamefile5_perm = [Outputdir,filesep,'Perm_ROI',sprintf('%05d',i),'_NetX2Y.nii'];
            end
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Res(i).res.ResultMap1;
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2y);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Res(i).res.ResultMap2;
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2x);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Res(i).res.ResultMap3;
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yT);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Res(i).res.ResultMap4;
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xT);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Res(i).res.ResultMap5;
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefile5);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Res(i).res.ResultMap_p1;
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2y_p);
            Dat = zeros(vinput(1).dim);
            Dat(find(DATMASK)) = Res(i).res.ResultMap_p2;
            DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2x_p);
            
            if PermMark
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).res.P_map1;
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2y_perm);
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).res.P_map2;
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2x_perm);
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).res.P_map3;
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yT_perm);
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).res.P_map4;
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xT_perm);
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).res.P_map5;
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefile5_perm);
            end
        end
    else
        for i = 1:size(ROIsignals,2)
            for i_ord = 1:RealCompPara.GCAorder
                
                outnamefilex2y = [Outputdir,filesep,'Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.nii'];
                outnamefiley2x = [Outputdir,filesep,'Coef_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.nii'];
                outnamefilex2yT = [Outputdir,filesep,'Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.nii'];
                outnamefiley2xT = [Outputdir,filesep,'Coef_T_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.nii'];
                outnamefilex2yZ = [Outputdir,filesep,'Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.nii'];
                outnamefiley2xZ = [Outputdir,filesep,'Coef_Z_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.nii'];
                outnamefilex2yP = [Outputdir,filesep,'Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.nii'];
                outnamefiley2xP = [Outputdir,filesep,'Coef_P_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.nii'];
                if PermMark                    
                    outnamefilex2yPerm = [Outputdir,filesep,'Coef_PermP_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_X2Y.nii'];
                    outnamefiley2xPerm = [Outputdir,filesep,'Coef_PermP_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'_Y2X.nii'];
                end
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).coef.ResultMap1{i_ord};
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2y);
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).coef.ResultMap2{i_ord};
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2x);
                
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).coef.res.T_T1{i_ord};
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yT);
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).coef.res.T_T2{i_ord};
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xT);
                
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).coef.res.Z_T1{i_ord};
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yZ);
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).coef.res.Z_T2{i_ord};
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xZ);
                
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).coef.res.P_T1{i_ord};
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yP);
                Dat = zeros(vinput(1).dim);
                Dat(find(DATMASK)) = Res(i).coef.res.P_T2{i_ord};
                DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xP);
                if PermMark
                    Dat = zeros(vinput(1).dim);
                    Dat(find(DATMASK)) = Res(i).coef.P_map1(i_ord,:);
                    DynamicBC_write_NIFTI(Dat,vinput(1),outnamefilex2yPerm);
                    Dat = zeros(vinput(1).dim);
                    Dat(find(DATMASK)) = Res(i).coef.P_map2(i_ord,:);
                    DynamicBC_write_NIFTI(Dat,vinput(1),outnamefiley2xPerm);
                end
            end
        end
    end
end
end


function indexs = mniroi(mnis,radius,Ttrans,dims,VOXELSIZE,DATAMASK)
coordinate = mni2cor(mnis, Ttrans);
sizes = round(radius./VOXELSIZE);
MATS = zeros(dims);
if coordinate(1)-sizes(1)>0
    Xrange(1) = coordinate(1)-sizes(1);
else
    Xrange(1) = 1;
end
if coordinate(1)+sizes(1)<=dims(1)
    Xrange(2) = coordinate(1)+sizes(1);
else
    Xrange(2) = dims(1);
end
%
if coordinate(2)-sizes(2)>0
    Yrange(1) = coordinate(2)-sizes(2);
else
    Yrange(1) = 1;
end
if coordinate(2)+sizes(2)<=dims(2)
    Yrange(2) = coordinate(2)+sizes(2);
else
    Yrange(2) = dims(2);
end
%
if coordinate(3)-sizes(3)>0
    Zrange(1) = coordinate(3)-sizes(3);
else
    Zrange(1) = 1;
end
if coordinate(3)+sizes(3)<=dims(3)
    Zrange(2) = coordinate(3)+sizes(3);
else
    Zrange(2) = dims(3);
end
MATS(Xrange(1):Xrange(2),Yrange(1):Yrange(2),Zrange(1):Zrange(2)) = 1;
MATS = MATS.*DATAMASK;
ind = find(MATS~=0);
[ix iy iz] = ind2sub(dims,ind);
exportnum = 1;
for i = 1:length(ind)
    if sqrt(((ix(i)-coordinate(1))*VOXELSIZE(1))^2+...
            ((iy(i)-coordinate(2))*VOXELSIZE(2))^2+...
            ((iz(i)-coordinate(3))*VOXELSIZE(3))^2)<radius
        indexs(exportnum) = ind(i);
        exportnum = exportnum+1;
    end
end
end