function BCCT_CON_computemain(Parameter)
pathfiles = which('BCCT.m');
RealCompPara.mod = 'volmap';
[path nam ext] = fileparts(pathfiles);
Outputdir = Parameter.Outputdir;
RealCompPara.Outputdir = Outputdir;
Inputdir = Parameter.Inputdir;
RealCompPara.Inputdir = Inputdir;
SeedROItype = Parameter.SeedROItype;
RealCompPara.SeedROItype = SeedROItype;
GlobalScale = Parameter.GlobalScal;
RealCompPara.GlobalScale = GlobalScale;

[vinput,datainput] = Dynamic_read_dir_NIFTI(Inputdir);
dims = vinput(1).dim;
Ttrans = vinput(1).mat;
VOXELSIZE = abs([vinput(1).mat(1,1),vinput(1).mat(2,2),vinput(1).mat(3,3)]);
RealCompPara.V = vinput(1);
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
    dynamicBC_Reslice(maskfiles,Maskimg0,vinput.dim,0,vinput.fname);
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
RealCompPara.DATMASK = DATMASK;
DATMASK = logical(DATMASK);
%% add in 2020-11-17
if GlobalScale % global mean scale
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
save([Outputdir,filesep,'RealCompPara.mat'],'RealCompPara')
for i = 1:size(ROIsignals,2)
    if i<10
        Outfilenametemp = fullfile(Outputdir,['R_ROI00000',num2str(i),'.nii']);
        OutfilenametempZ = fullfile(Outputdir,['Z_ROI00000',num2str(i),'.nii']);
        OutfilenametempP = fullfile(Outputdir,['P_ROI00000',num2str(i),'.nii']);
    elseif i<100
        Outfilenametemp = fullfile(Outputdir,['R_ROI0000',num2str(i),'.nii']);
        OutfilenametempZ = fullfile(Outputdir,['Z_ROI0000',num2str(i),'.nii']);
        OutfilenametempP = fullfile(Outputdir,['P_ROI0000',num2str(i),'.nii']);
    elseif i<1000
        Outfilenametemp = fullfile(Outputdir,['R_ROI000',num2str(i),'.nii']);
        OutfilenametempZ = fullfile(Outputdir,['Z_ROI000',num2str(i),'.nii']);
        OutfilenametempP = fullfile(Outputdir,['P_ROI000',num2str(i),'.nii']);
    else
        Outfilenametemp = fullfile(Outputdir,['R_ROI00',num2str(i),'.nii']);
        OutfilenametempZ = fullfile(Outputdir,['Z_ROI00',num2str(i),'.nii']);
        OutfilenametempP = fullfile(Outputdir,['P_ROI00',num2str(i),'.nii']);
    end
    if COVcond==1
        [R P] = partialcorr(maskedSignal,ROIsignals(:,i),COV);
        DF_E = size(ROIsignals,1)-2-size(COV,2);
        [Z P2] = AS_TFRtoZ(R,'R',DF_E,[]);
    else
        [R P] = corr(maskedSignal,ROIsignals(:,i));
        DF_E = size(ROIsignals,1)-2;
        [Z P2] = AS_TFRtoZ(R,'R',DF_E,[]);
    end
    Rmaps = zeros(dims);
    Rmaps(find(DATMASK)) = R;
    vmasknew = vmask;
    vmasknew.descrip=sprintf('{R_[%.1f]}',DF_E);
    DynamicBC_write_NIFTI(Rmaps,vmask,Outfilenametemp);
    Zmaps = zeros(dims);
    Zmaps(find(DATMASK)) = Z;
    DynamicBC_write_NIFTI(Zmaps,vmask,OutfilenametempZ);
    Pmaps = zeros(dims);
    Pmaps(find(DATMASK)) = P;
    DynamicBC_write_NIFTI(Pmaps,vmask,OutfilenametempP);
end
%% This part will consider in master version.
% if strcmp(Parameter.cal_types,'Yes')
%     for i = 1:size(ROIsignals,2)
%         if i<10
%             Outfilenametemp = fullfile(Outputdir,['Scatter_ROI00000',num2str(i),'.nii']);
%             OutfilenametempP = fullfile(Outputdir,['Scatter_P_ROI00000',num2str(i),'.nii']);
%             Outfilemat = fullfile(Outputdir,['Scatter_ROI00000',num2str(i),'.mat']);
%         elseif i<100
%             Outfilenametemp = fullfile(Outputdir,['Scatter_ROI0000',num2str(i),'.nii']);
%             OutfilenametempP = fullfile(Outputdir,['Scatter_P_ROI0000',num2str(i),'.nii']);
%             Outfilemat = fullfile(Outputdir,['Scatter_ROI0000',num2str(i),'.mat']);
%         elseif i<1000
%             Outfilenametemp = fullfile(Outputdir,['Scatter_ROI000',num2str(i),'.nii']);
%             OutfilenametempP = fullfile(Outputdir,['Scatter_P_ROI000',num2str(i),'.nii']);
%             Outfilemat = fullfile(Outputdir,['Scatter_PROI000',num2str(i),'.mat']);
%         else
%             Outfilenametemp = fullfile(Outputdir,['Scatter_ROI00',num2str(i),'.nii']);
%             OutfilenametempP = fullfile(Outputdir,['Scatter_P_ROI00',num2str(i),'.nii']);
%             Outfilemat = fullfile(Outputdir,['Scatter_ROI00',num2str(i),'.mat']);
%         end
% %         if COVcond==1
% %             [R P] = partialcorr(maskedSignal,ROIsignals(:,i),COV);
% %         else
% %             [R P] = corr(maskedSignal,ROIsignals(:,i));
% %         end
%         [R P labs] = scatcorr(maskedSignal,ROIsignals(:,i));
%         Rmaps = zeros(dims);
%         Rmaps(find(DATMASK)) = R;
%         DynamicBC_write_NIFTI(Rmaps,vmask,Outfilenametemp);
%         Pmaps = zeros(dims);
%         Pmaps(find(DATMASK)) = P;
%         DynamicBC_write_NIFTI(Pmaps,vmask,OutfilenametempP);
%         save(Outfilemat,'labs');
%     end
% end
%%
disp('Compute Finished!')
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