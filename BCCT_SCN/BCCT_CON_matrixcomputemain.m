function BCCT_CON_matrixcomputemain(Parameter)
pathfiles = which('BCCT.m');
RealCompPara.mod = 'volroi';
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


if SeedROItype<3
    [vinput,datainput] = Dynamic_read_dir_NIFTI(Inputdir);
    dims = vinput.dim;
    Ttrans = vinput.mat;
    VOXELSIZE = abs([vinput.mat(1,1),vinput.mat(2,2),vinput.mat(3,3)]);
    if strcmp(Parameter.masks,'Defaults')
        maskfiles = fullfile(path,'template','grey.nii');
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
save(fullfile(Outputdir,'RealCompPara.mat'),'RealCompPara');
%% Consider in master version.
% cal_types = Parameter.cal_types;
% if strcmp(cal_types,'Yes')
%     R = zeros(nrois);
%     P = ones(nrois);
%     kind = 1;
%     for i = 1:nrois-1
%         for j = i+1:nrois
%             Rsig1 = ROIsignals(:,i);
%             Rsig2 = ROIsignals(:,j);
%             [Pi p ol] = Shepherd(Rsig1,Rsig2);
%             Rsigtemp = [Rsig1,Rsig2];
%             Rsigtemp(find(ol),:) = [];
%             [rtemp ptemp] = corr(Rsigtemp); 
%             R(i,j) = rtemp(2);
%             P(i,j) = ptemp(2);
%             R(j,i) = rtemp(2);
%             P(j,i) = ptemp(2);
%             Sellab(:,kind) = sparse(ol);
%             Selind{kind} = [i,j];
%             kind = kind+1;
%         end
%     end
%     Outputfile = fullfile(Outputdir,'Scatter_R_Pres.mat');
%     save(Outputfile,'R','P','Sellab','Selind','nrois');
% end
%%
if COVcond==1
    if Partinfo
        R = zeros(size(ROIsignals,2));
        P = ones(size(ROIsignals,2));
        for i = 1:size(ROIsignals,2)-1
            for j = i+1:size(ROIsignals,2)
                sigind = 1:size(ROIsignals,2);
                sigind([i,j]) = [];
                COVnew = [COV,ROIsignals(:,sigind)];
                [r p] = partialcorr(ROIsignals(:,i),ROIsignals(:,j),COVnew);
                R(i,j) = r;
                P(i,j) = p;
                R(j,i) = r;
                P(j,i) = p;
            end
        end
        DF_E = size(ROIsignals,1)-2-size(COVnew,2);
        [Z,PZ] = AS_TFRtoZ(R,'R',DF_E,[]); 
    else
        [R P] = partialcorr(ROIsignals,COV);
        DF_E = size(ROIsignals,1)-2-size(COV,2);
        [Z, PZ] = AS_TFRtoZ(R,'R',DF_E,[]);        
    end
else
    if Partinfo        
        R = zeros(size(ROIsignals,2));
        P = ones(size(ROIsignals,2));
        for i = 1:size(ROIsignals,2)-1
            for j = i+1:size(ROIsignals,2)
                sigind = 1:size(ROIsignals,2);
                sigind([i,j]) = [];
                COVnew = ROIsignals(:,sigind);
                [r p] = partialcorr(ROIsignals(:,i),ROIsignals(:,j),COVnew);
                R(i,j) = r;
                P(i,j) = p;
                R(j,i) = r;
                P(j,i) = p;
            end
        end
        DF_E = size(ROIsignals,1)-2-size(COVnew,2);
        [Z,PZ] = AS_TFRtoZ(R,'R',DF_E,[]); 
    else
        [R P] = corr(ROIsignals);
        DF_E = size(ROIsignals,1)-2;
        [Z, PZ] = AS_TFRtoZ(R,'R',DF_E,[]);        
    end
end
Outputfile = fullfile(Outputdir,'R_Pres.mat');
save(Outputfile,'R','P','Z', 'PZ','DF_E')
disp('Compute Finished!')
end