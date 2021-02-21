function BCCT_WTA_Stat_main_ver2(Parameter)
% Face a problem for the number of connects
% In later version, we would like to deal with it for the small ROM
% computer.

% Change this value to speed up when you have large ROM
Nslice = 1000;
%

Outdir = Parameter.Outdir;
Indir1 = Parameter.Input1;
Indir2 = Parameter.Input2;
Incalmat1 = load(fullfile(Indir1,'SetUpparameter.mat')); 
Incalmat2 = load(fullfile(Indir2,'SetUpparameter.mat'));
COVcond1 = Incalmat1.Parameter.covs; % g1 协变量情况
COVcond2 = Incalmat2.Parameter.covs; % g2 协变量情况
Calmethod1 = Incalmat1.Parameter.methodused; % G1 计算方法 1 pearson 2 partialcorr
Calmethod2 = Incalmat2.Parameter.methodused; % G2 计算方法 1 pearson 2 partialcorr 
if COVcond1~=COVcond2
    error('The two groups used different cov conditions');
end
if Calmethod1~=Calmethod2
    error('The tow groups used different calculate methods');
end
COVCOND = COVcond1;
PermMethod = Parameter.permtype; % 统计方法： 1 permutation 2 interaction 

computeval1 = load(fullfile(Indir1,'computeval.mat'));
computeval2 = load(fullfile(Indir2,'computeval.mat'));
OrigR1 = computeval1.r; % G1 no scat的r值
OrigR2 = computeval2.r; % G2 no scat的r值
DELTR = OrigR1-OrigR2; % 差异值，主要用于permutation

[OrigZ1,~] = AS_TFRtoZ(computeval1.r,'R',computeval1.DF_E,[]);
[OrigZ2,~] = AS_TFRtoZ(computeval2.r,'R',computeval2.DF_E,[]);
DELTZ = OrigZ1-OrigZ2;
DF_E1 = computeval1.DF_E;
DF_E2 = computeval2.DF_E;

SIG1 = load(fullfile(Indir1,'Signals.mat')); % 读取原始信号
SIG2 = load(fullfile(Indir2,'Signals.mat')); % 读取原始信号
datTarget1 = SIG1.datTarget; % G1 target信号（原始）
datTarget2 = SIG2.datTarget; % G2 target信号（原始）
datSeed1 = SIG1.datSeed; % G1 seed信号（原始）
datSeed2 = SIG2.datSeed; % G2 seed信号（原始）
NROI = size(datSeed1,1);

[vtargROI dtargROI] = Dynamic_read_dir_NIFTI(Incalmat1.Parameter.TargetROI);
valtargROI = unique(dtargROI);
MulTarROI = 0;
if length(valtargROI)>2 % multidefinedROI
    MulTarROI = 1;
end

if PermMethod==1 % permutation 
    Permnum = Parameter.permnum; % permutation数目
    N1 = size(datTarget1,2); % G1 被试数（时间点长度）
    N2 = size(datTarget2,2); % G2 被试数（时间点长度）
    Ntotal = N1+N2; % 总体被试数目（时间点长度）
    dattarget = [datTarget1,datTarget2]'; % 拼接的target数据
    if COVCOND==0 % no cov 
%% cond 1
        if Calmethod1==1 % pearson & nocov
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]'; % 第i个ROI的seed信号拼接
%                 deltr = zeros(Permnum,size(datTarget1,1)); % 预设值delta-p值矩阵
%                 Rtemp1 = zeros(Permnum,size(datTarget1,1));  
%                 Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                % for the small ROI computer, slice comp will help it.
%                 Nslice = 1000;
                Totalnum = size(datTarget1,1);
                Numcal = ceil(Totalnum/Nslice);
                for iN = 1:Numcal
                    disp(['round ',num2str(iN),', total ',num2str(Numcal)]);
                    tic
                    if iN<Numcal
                        Sliceind = 1+Nslice*(iN-1):Nslice*iN;
                    else
                        Sliceind = 1+Nslice*(iN-1):Totalnum;
                    end
                    parfor iperm = 1:Permnum
                        Nrandp = randperm(Ntotal);
                        seedsig1 = sigseedt(Nrandp(1:N1));
                        seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                        targetsig1 = dattarget(Nrandp(1:N1),Sliceind);  % 前N1个target
                        targetsig2 = dattarget(Nrandp(N1+1:Ntotal),Sliceind); % 后N2个target
                        [r1 p1] = corr(seedsig1,targetsig1); %相关分析
                        [r2 p2] = corr(seedsig2,targetsig2); %相关分析
                        deltr(iperm,:) = r1-r2; % delta-p值录入
                        Rtemp1(iperm,:) = r1;
                        Rtemp2(iperm,:) = r2;
                    end
                    toc
                    Rtemp1(isnan(Rtemp1)) = 0;
                    Rtemp2(isnan(Rtemp2)) = 0;
                    [Ztemp1,ZPtemp1] = AS_TFRtoZ(Rtemp1,'R',DF_E1,[]);
                    [Ztemp2,ZPtemp2] = AS_TFRtoZ(Rtemp2,'R',DF_E2,[]);
                    save([Outdir,filesep,'ExtRval_',num2str(i),'_round',num2str(iN),'(total)',num2str(Numcal),'.mat'],'Rtemp1','Rtemp2','Sliceind','Ztemp1','Ztemp2');
                    [mu,sig,~,sigci] = normfit(deltr); %正态分布
                    P = normcdf(DELTR(i,Sliceind),mu,sig);
                    Pval(i,Sliceind) = P;
                    [mu,sig,~,sigci] = normfit(Ztemp1-Ztemp2);
                    ZP = normcdf(DELTZ(i,Sliceind),mu,sig);
                    ZPval(i,Sliceind) = ZP;
                    clear Rtemp1 Rtemp2 Ztemp1 Ztemp2 Sliceind mu sig sigci deltr
                    toc
                end
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
                ZPtemp = zeros(vmat.dim);
                ZPtemp(INDEXS) = ZPval(i,:);
                DynamicBC_write_NIFTI(ZPtemp,vmat,OutPnameZ);
            end
%% Condition 2
        else % particalcorr & nocov
            for i = 1:NROI % 逐个ROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]';  % 拼接种子点
%                 deltr = zeros(Permnum,size(datTarget1,1));  % 定义r差异矩阵
%                 Rtemp1 = zeros(Permnum,size(datTarget1,1));  
%                 Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                INDUSED = 1:NROI;INDUSED(i) = []; % 定义协变量
                covt = [datSeed1(INDUSED,:),datSeed2(INDUSED,:)]'; % 协变量拼接
%                 Nslice = 1000;
                Totalnum = size(datTarget1,1);
                Numcal = ceil(Totalnum/Nslice);
                for iN = 1:Numcal
                    disp(['round ',num2str(iN),', total ',num2str(Numcal)]);
                    tic
                    if iN<Numcal
                        Sliceind = 1+Nslice*(iN-1):Nslice*iN;
                    else
                        Sliceind = 1+Nslice*(iN-1):Totalnum;
                    end
                    parfor iperm = 1:Permnum
                        Nrandp = randperm(Ntotal); % 随机排序
                        seedsig1 = sigseedt(Nrandp(1:N1)); % 前N1个Seed
                        seedsig2 = sigseedt(Nrandp(N1+1:Ntotal)); % 后N2个seed
                        targetsig1 = dattarget(Nrandp(1:N1),Sliceind); % 前N1个target
                        targetsig2 = dattarget(Nrandp(N1+1:Ntotal),Sliceind); % 后N2个target
                        covt1 = covt(Nrandp(1:N1),:); % 前N1个COV
                        covt2 = covt(Nrandp(N1+1:Ntotal),:); % 后N2个COV
                        [r1 p1] = partialcorr(seedsig1,targetsig1,covt1); % 偏相关
                        [r2 p2] = partialcorr(seedsig2,targetsig2,covt2); % 偏相关
                        deltr(iperm,:) = r1-r2; % r值差异录入
                        Rtemp1(iperm,:) = r1;
                        Rtemp2(iperm,:) = r2;
                    end
                    Rtemp1(isnan(Rtemp1)) = 0;
                    Rtemp2(isnan(Rtemp2)) = 0;
                    [Ztemp1,ZPtemp1] = AS_TFRtoZ(Rtemp1,'R',DF_E1,[]);
                    [Ztemp2,ZPtemp2] = AS_TFRtoZ(Rtemp2,'R',DF_E2,[]);
                    save([Outdir,filesep,'ExtRval_',num2str(i),'_round',num2str(iN),'(total)',num2str(Numcal),'.mat'],'Rtemp1','Rtemp2','Sliceind','Ztemp1','Ztemp2');
                    [mu,sig,~,sigci] = normfit(deltr); %正态分布
                    P = normcdf(DELTR(i,Sliceind),mu,sig);
                    Pval(i,Sliceind) = P;
                    [mu,sig,~,sigci] = normfit(Ztemp1-Ztemp2);
                    ZP = normcdf(DELTZ(i,Sliceind),mu,sig);
                    ZPval(i,Sliceind) = ZP;
                    clear Rtemp1 Rtemp2 Ztemp1 Ztemp2 Sliceind mu sig sigci deltr
                    toc
                end
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
                ZPtemp = zeros(vmat.dim);
                ZPtemp(INDEXS) = ZPval(i,:);
                DynamicBC_write_NIFTI(ZPtemp,vmat,OutPnameZ);
            end
        end
    else % with cov
        if Calmethod1==1 % pearson & withcov
%% condition 3
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]'; % 第i个ROI
%                 deltr = zeros(Permnum,size(datTarget1,1)); % 定义delta-r
%                 Rtemp1 = zeros(Permnum,size(datTarget1,1));  
%                 Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                covt = [SIG1.COV;SIG2.COV]; % 协变量
                
%                 Nslice = 1000;
                Totalnum = size(datTarget1,1);
                Numcal = ceil(Totalnum/Nslice);
                for iN = 1:Numcal
                    disp(['round ',num2str(iN),', total ',num2str(Numcal)]);
                    tic
                    if iN<Numcal
                        Sliceind = 1+Nslice*(iN-1):Nslice*iN;
                    else
                        Sliceind = 1+Nslice*(iN-1):Totalnum;
                    end
                    parfor iperm = 1:Permnum
                        Nrandp = randperm(Ntotal);
                        seedsig1 = sigseedt(Nrandp(1:N1));
                        seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                        targetsig1 = dattarget(Nrandp(1:N1),Sliceind);
                        targetsig2 = dattarget(Nrandp(N1+1:Ntotal),Sliceind);
                        covt1 = covt(Nrandp(1:N1),:); % 前N1个COV
                        covt2 = covt(Nrandp(N1+1:Ntotal),:); % 后N2个COV
                        [r1 p1] = partialcorr(seedsig1,targetsig1,covt1); % 偏相关
                        [r2 p2] = partialcorr(seedsig2,targetsig2,covt2); % 偏相关
                        deltr(iperm,:) = r1-r2; % r值差异录入
                        Rtemp1(iperm,:) = r1;
                        Rtemp2(iperm,:) = r2;
                    end
                    Rtemp1(isnan(Rtemp1)) = 0;
                    Rtemp2(isnan(Rtemp2)) = 0;
                    [Ztemp1,ZPtemp1] = AS_TFRtoZ(Rtemp1,'R',DF_E1,[]);
                    [Ztemp2,ZPtemp2] = AS_TFRtoZ(Rtemp2,'R',DF_E2,[]);
                    save([Outdir,filesep,'ExtRval_',num2str(i),'_round',num2str(iN),'(total)',num2str(Numcal),'.mat'],'Rtemp1','Rtemp2','Sliceind','Ztemp1','Ztemp2');
                    [mu,sig,~,sigci] = normfit(deltr); %正态分布
                    P = normcdf(DELTR(i,Sliceind),mu,sig);
                    Pval(i,Sliceind) = P;
                    [mu,sig,~,sigci] = normfit(Ztemp1-Ztemp2);
                    ZP = normcdf(DELTZ(i,Sliceind),mu,sig);
                    ZPval(i,Sliceind) = ZP;
                    clear Rtemp1 Rtemp2 Ztemp1 Ztemp2 Sliceind mu sig sigci deltr
                    toc
                end
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
                ZPtemp = zeros(vmat.dim);
                ZPtemp(INDEXS) = ZPval(i,:);
                DynamicBC_write_NIFTI(ZPtemp,vmat,OutPnameZ);
            end
        else % partialcorr with cov
%% condition 4
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]';
%                 deltr = zeros(Permnum,size(datTarget1,1));
%                 Rtemp1 = zeros(Permnum,size(datTarget1,1));  
%                 Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                covtCOV = [SIG1.COV;SIG2.COV];
                INDUSED = 1:NROI;INDUSED(i) = [];
                covt = [covtCOV,[datSeed1(INDUSED,:),datSeed2(INDUSED,:)]']; % cov由COV和其他seed信号构成
                
%                 Nslice = 1000;
                Totalnum = size(datTarget1,1);
                Numcal = ceil(Totalnum/Nslice);
                for iN = 1:Numcal
                    disp(['round ',num2str(iN),', total ',num2str(Numcal)]);
                    tic
                    if iN<Numcal
                        Sliceind = 1+Nslice*(iN-1):Nslice*iN;
                    else
                        Sliceind = 1+Nslice*(iN-1):Totalnum;
                    end
                    parfor iperm = 1:Permnum
                        Nrandp = randperm(Ntotal);
                        seedsig1 = sigseedt(Nrandp(1:N1));
                        seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                        targetsig1 = dattarget(Nrandp(1:N1),Sliceind);
                        targetsig2 = dattarget(Nrandp(N1+1:Ntotal),Sliceind);
                        covt1 = covt(Nrandp(1:N1),:); % 前N1个COV
                        covt2 = covt(Nrandp(N1+1:Ntotal),:); % 后N2个COV
                        [r1 p1] = partialcorr(seedsig1,targetsig1,covt1); % 偏相关
                        [r2 p2] = partialcorr(seedsig2,targetsig2,covt2); % 偏相关
                        deltr(iperm,:) = r1-r2; % r值差异录入
                        Rtemp1(iperm,:) = r1;
                        Rtemp2(iperm,:) = r2;
                    end
                    Rtemp1(isnan(Rtemp1)) = 0;
                    Rtemp2(isnan(Rtemp2)) = 0;
                    [Ztemp1,ZPtemp1] = AS_TFRtoZ(Rtemp1,'R',DF_E1,[]);
                    [Ztemp2,ZPtemp2] = AS_TFRtoZ(Rtemp2,'R',DF_E2,[]);
                    save([Outdir,filesep,'ExtRval_',num2str(i),'_round',num2str(iN),'(total)',num2str(Numcal),'.mat'],'Rtemp1','Rtemp2','Sliceind','Ztemp1','Ztemp2');
                    [mu,sig,~,sigci] = normfit(deltr); %正态分布
                    P = normcdf(DELTR(i,Sliceind),mu,sig);
                    Pval(i,Sliceind) = P;
                    [mu,sig,~,sigci] = normfit(Ztemp1-Ztemp2);
                    ZP = normcdf(DELTZ(i,Sliceind),mu,sig);
                    ZPval(i,Sliceind) = ZP;
                    clear Rtemp1 Rtemp2 Ztemp1 Ztemp2 Sliceind mu sig sigci deltr
                    toc
                end   
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                    OutPnameZ = fullfile(Outdir,['ZPval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
                ZPtemp = zeros(vmat.dim);
                ZPtemp(INDEXS) = ZPval(i,:);
                DynamicBC_write_NIFTI(ZPtemp,vmat,OutPnameZ);
            end
        end
    end
else % interaction
    
end
%% This part is only for the permutation, test the number of voxel significant
% this is  for pos
if PermMethod==1 % permutation    
    zn1 = zeros(Permnum,Totalnum);
    zn2 = zeros(Permnum,Totalnum);
    n1 = zeros(Permnum,Totalnum);
    n2 = zeros(Permnum,Totalnum);
    for i = 1:Numcal
        for j = 1:NROI
            %             ,'Rtemp1','Rtemp2','Sliceind','Ztemp1','Ztemp2'
            clear Rtemp1 Rtemp2 Sliceind Ztemp1 Ztemp2 Rg1 Rg2
            Rg1 = load([Outdir,filesep,'ExtRval_',num2str(j),'_round',num2str(i),'(total)',num2str(Numcal),'.mat']);
            Rg2 = load([Outdir,filesep,'ExtRval_',num2str(j),'_round',num2str(i),'(total)',num2str(Numcal),'.mat']);
            RvalComp1(:,:,j) = Rg1.Rtemp1;
            RvalComp2(:,:,j) = Rg2.Rtemp2;
            ZvalComp1(:,:,j) = Rg1.Ztemp1;
            ZvalComp2(:,:,j) = Rg2.Ztemp2;
            Sliceind = Rg1.Sliceind;
        end
        clear N1 N2 ZN1 ZN2 
        for k = 1:length(Sliceind)
            RvalComp1temp = squeeze(RvalComp1(:,k,:));
            RvalComp2temp = squeeze(RvalComp2(:,k,:));
            [maxv, imaxv] = max(RvalComp1temp');
            N1(:,k) = imaxv;
            [maxv, imaxv] = max(RvalComp2temp');
            N2(:,k) = imaxv;
            ZvalComp1temp = squeeze(ZvalComp1(:,k,:));
            ZvalComp2temp = squeeze(ZvalComp2(:,k,:));
            [maxv, imaxv] = max(ZvalComp1temp');
            ZN1(:,k) = imaxv;
            [maxv, imaxv] = max(ZvalComp2temp');
            ZN2(:,k) = imaxv;
        end
        n1(:,Sliceind) = N1;
        n2(:,Sliceind) = N2;
        zn1(:,Sliceind) = ZN1;
        zn2(:,Sliceind) = ZN2;
        clear N1 N2 ZN1 ZN2 RvalComp1 RvalComp2 ZvalComp1 ZvalComp2
    end
    for i = 1:Permnum
        for j = 1:NROI
            ExtN1(i,j) = length(find(n1(i,:)==j));
            ExtN2(i,j) = length(find(n2(i,:)==j));
            ExtZN1(i,j) = length(find(zn1(i,:)==j));
            ExtZN2(i,j) = length(find(zn2(i,:)==j));
        end
    end
    save([Outdir,filesep,'ConnectNum.mat'],'ExtN1','ExtN2','ExtZN1','ExtZN2')
    if MulTarROI==1
        DTARGROI = dtargROI(INDEXS);
        for iroi = 1:length(valtargROI)-1
            indtempseproi = find(DTARGROI==valtargROI(iroi+1));
            for i = 1:Permnum
                for j = 1:NROI
                    ExtN1ROI(i,j,iroi) = length(find(n1(i,indtempseproi)==j));
                    ExtN2ROI(i,j,iroi) = length(find(n2(i,indtempseproi)==j));
                    ExtZN1ROI(i,j,iroi) = length(find(zn1(i,indtempseproi)==j));
                    ExtZN2ROI(i,j,iroi) = length(find(zn2(i,indtempseproi)==j));
                end
            end
        end
        save([Outdir,filesep,'ConnectNumROIS.mat'],'ExtN1ROI','ExtN2ROI','ExtZN1ROI','ExtZN2ROI')
    end
end
% abs
if PermMethod==1 % permutation    
    zn1 = zeros(Permnum,Totalnum);
    zn2 = zeros(Permnum,Totalnum);
    n1 = zeros(Permnum,Totalnum);
    n2 = zeros(Permnum,Totalnum);
    for i = 1:Numcal
        for j = 1:NROI
            %             ,'Rtemp1','Rtemp2','Sliceind','Ztemp1','Ztemp2'
            clear Rtemp1 Rtemp2 Sliceind Ztemp1 Ztemp2 Rg1 Rg2
            Rg1 = load([Outdir,filesep,'ExtRval_',num2str(j),'_round',num2str(i),'(total)',num2str(Numcal),'.mat']);
            Rg2 = load([Outdir,filesep,'ExtRval_',num2str(j),'_round',num2str(i),'(total)',num2str(Numcal),'.mat']);
            RvalComp1(:,:,j) = abs(Rg1.Rtemp1);
            RvalComp2(:,:,j) = abs(Rg2.Rtemp2);
            ZvalComp1(:,:,j) = abs(Rg1.Ztemp1);
            ZvalComp2(:,:,j) = abs(Rg2.Ztemp2);
            Sliceind = Rg1.Sliceind;
        end
        clear N1 N2 ZN1 ZN2 
        for k = 1:length(Sliceind)
            RvalComp1temp = squeeze(RvalComp1(:,k,:));
            RvalComp2temp = squeeze(RvalComp2(:,k,:));
            [maxv, imaxv] = max(RvalComp1temp');
            N1(:,k) = imaxv;
            [maxv imaxv] = max(RvalComp2temp');
            N2(:,k) = imaxv;
            ZvalComp1temp = squeeze(ZvalComp1(:,k,:));
            ZvalComp2temp = squeeze(ZvalComp2(:,k,:));
            [maxv imaxv] = max(ZvalComp1temp');
            ZN1(:,k) = imaxv;
            [maxv imaxv] = max(ZvalComp2temp');
            ZN2(:,k) = imaxv;
        end
        n1(:,Sliceind) = N1;
        n2(:,Sliceind) = N2;
        zn1(:,Sliceind) = ZN1;
        zn2(:,Sliceind) = ZN2;
        clear N1 N2 ZN1 ZN2 RvalComp1 RvalComp2 ZvalComp1 ZvalComp2
    end
    for i = 1:Permnum
        for j = 1:NROI
            ExtN1(i,j) = length(find(n1(i,:)==j));
            ExtN2(i,j) = length(find(n2(i,:)==j));
            ExtZN1(i,j) = length(find(zn1(i,:)==j));
            ExtZN2(i,j) = length(find(zn2(i,:)==j));
        end
    end
    save([Outdir,filesep,'ABS_ConnectNum.mat'],'ExtN1','ExtN2','ExtZN1','ExtZN2')
    if MulTarROI==1
        DTARGROI = dtargROI(INDEXS);
        for iroi = 1:length(valtargROI)-1
            indtempseproi = find(DTARGROI==valtargROI(iroi+1));
            for i = 1:Permnum
                for j = 1:NROI
                    ExtN1ROI(i,j,iroi) = length(find(n1(i,indtempseproi)==j));
                    ExtN2ROI(i,j,iroi) = length(find(n2(i,indtempseproi)==j));
                    ExtZN1ROI(i,j,iroi) = length(find(zn1(i,indtempseproi)==j));
                    ExtZN2ROI(i,j,iroi) = length(find(zn2(i,indtempseproi)==j));
                end
            end
        end
        save([Outdir,filesep,'ABS_ConnectNumROIS.mat'],'ExtN1ROI','ExtN2ROI','ExtZN1ROI','ExtZN2ROI')
    end
end

end
% [vtargROI dtargROI] = Dynamic_read_dir_NIFTI(Incalmat1.Parameter.TargetROI);
% valtargROI = unique(dtargROI);
% MulTarROI = 0;
% if length(valtargROI)>2 % multidefinedROI
%     MulTarROI = 1;
% end