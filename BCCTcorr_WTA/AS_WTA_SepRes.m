function AS_WTA_SepRes
indir = uigetdir(pwd,'OrigData');
outdir = uigetdir(pwd,'Outputdir');
[TargetROI,FilPa,Filext] = uigetfile({'*.nii';'*.img'},'TargetROI');
[vtar,dtar] = Dynamic_read_dir_NIFTI(fullfile(FilPa,TargetROI));
NIINAMELIST = dir([indir,filesep,'*.nii']);
dtarind = unique(dtar);
for i = 1:length(NIINAMELIST)
    [VT,DT] = Dynamic_read_dir_NIFTI(fullfile(indir,NIINAMELIST(i).name));
    for j = 1:length(dtarind)-1
        DATO = DT.*(dtar==dtarind(j+1));
        DATT = reshape(DATO,VT.dim(1),VT.dim(2),VT.dim(3));
        DynamicBC_write_NIFTI(DATT,VT,fullfile(outdir,['ROI_',num2str(j),'_',NIINAMELIST(i).name]));
    end
end
