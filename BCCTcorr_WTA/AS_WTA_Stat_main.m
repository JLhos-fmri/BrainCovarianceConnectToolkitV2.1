function AS_WTA_Stat_main(Parameter)
% Face a problem for the number of connects
% In later version, we would like to deal with it for the small ROM
% computer.
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
SIG1 = load(fullfile(Indir1,'Signals.mat')); % 读取原始信号
SIG2 = load(fullfile(Indir2,'Signals.mat')); % 读取原始信号
datTarget1 = SIG1.datTarget; % G1 target信号（原始）
datTarget2 = SIG2.datTarget; % G2 target信号（原始）
datSeed1 = SIG1.datSeed; % G1 seed信号（原始）
datSeed2 = SIG2.datSeed; % G2 seed信号（原始）
NROI = size(datSeed1,1);
if PermMethod==1 % permutation 
    Permnum = Parameter.permnum; % permutation数目
    N1 = size(datTarget1,2); % G1 被试数（时间点长度）
    N2 = size(datTarget2,2); % G2 被试数（时间点长度）
    Ntotal = N1+N2; % 总体被试数目（时间点长度）
    dattarget = [datTarget1,datTarget2]'; % 拼接的target数据
    if COVCOND==0 % no cov 
        if Calmethod1==1 % pearson & nocov
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]'; % 第i个ROI的seed信号拼接
                deltr = zeros(Permnum,size(datTarget1,1)); % 预设值delta-p值矩阵
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                parfor iperm = 1:Permnum % 逐次置换
                    Nrandp = randperm(Ntotal); % 重排
                    seedsig1 = sigseedt(Nrandp(1:N1));  % 前N1个Seed
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal)); % 后N2个Seed
                    targetsig1 = dattarget(Nrandp(1:N1),:);  % 前N1个target
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:); % 后N2个target
                    [r1 p1] = corr(seedsig1,targetsig1); %相关分析
                    [r2 p2] = corr(seedsig2,targetsig2); %相关分析
                    deltr(iperm,:) = r1-r2; % delta-p值录入
                    Rtemp1(iperm,:) = r1; 
                    Rtemp2(iperm,:) = r2;
                end
                save([Outdir,filesep,'ExtRval_',num2str(i),'.mat'],'Rtemp1','Rtemp2');
                [mu,sig,~,sigci] = normfit(deltr); %正态分布
                P = normcdf(DELTR(i,:),mu,sig); % 计算p值（<0.05为G1<G2，>0.95为G1>G2)
                Pval(i,:) = P; % 录入p值
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        else % particalcorr & nocov
            for i = 1:NROI % 逐个ROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]';  % 拼接种子点
                deltr = zeros(Permnum,size(datTarget1,1));  % 定义r差异矩阵
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                INDUSED = 1:NROI;INDUSED(i) = []; % 定义协变量
                covt = [datSeed1(INDUSED,:),datSeed2(INDUSED,:)]'; % 协变量拼接
                parfor iperm = 1:Permnum % 逐次置换
                    Nrandp = randperm(Ntotal); % 随机排序
                    seedsig1 = sigseedt(Nrandp(1:N1)); % 前N1个Seed
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal)); % 后N2个seed
                    targetsig1 = dattarget(Nrandp(1:N1),:); % 前N1个target
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:); % 后N2个target
                    covt1 = covt(Nrandp(1:N1),:); % 前N1个COV
                    covt2 = covt(Nrandp(N1+1:Ntotal),:); % 后N2个COV
                    [r1 p1] = partialcorr(seedsig1,targetsig1,covt1); % 偏相关
                    [r2 p2] = partialcorr(seedsig2,targetsig2,covt2); % 偏相关
                    deltr(iperm,:) = r1-r2; % r值差异录入
                    Rtemp1(iperm,:) = r1; 
                    Rtemp2(iperm,:) = r2;
                end
                save([Outdir,filesep,'ExtRval_',num2str(i),'.mat'],'Rtemp1','Rtemp2');
                [mu,sig,~,sigci] = normfit(deltr); % 正态分布
                P = normcdf(DELTR(i,:),mu,sig); % p值计算
                Pval(i,:) = P; %p值录入
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        end
    else % with cov
        if Calmethod1==1 % pearson & withcov
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]'; % 第i个ROI
                deltr = zeros(Permnum,size(datTarget1,1)); % 定义delta-r
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                covt = [SIG1.COV;SIG2.COV]; % 协变量
                parfor iperm = 1:Permnum % 逐次置换
                    Nrandp = randperm(Ntotal);
                    seedsig1 = sigseedt(Nrandp(1:N1));
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                    targetsig1 = dattarget(Nrandp(1:N1),:);
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:);
                    covt1 = covt(Nrandp(1:N1),:);
                    covt2 = covt(Nrandp(N1+1:Ntotal),:);
                    [r1 p1] = partialcorr(seedsig1,targetsig1,covt1);
                    [r2 p2] = partialcorr(seedsig2,targetsig2,covt2);
                    deltr(iperm,:) = r1-r2;
                    
                    Rtemp1(iperm,:) = r1; 
                    Rtemp2(iperm,:) = r2;
                end
                save([Outdir,filesep,'ExtRval_',num2str(i),'.mat'],'Rtemp1','Rtemp2');
                [mu,sig,~,sigci] = normfit(deltr);
                P = normcdf(DELTR(i,:),mu,sig);
                Pval(i,:) = P;
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        else % partialcorr with cov
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]';
                deltr = zeros(Permnum,size(datTarget1,1));
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                covtCOV = [SIG1.COV;SIG2.COV];
                INDUSED = 1:NROI;INDUSED(i) = [];
                covt = [covtCOV,[datSeed1(INDUSED,:),datSeed2(INDUSED,:)]']; % cov由COV和其他seed信号构成
                parfor iperm = 1:Permnum
                    Nrandp = randperm(Ntotal);
                    seedsig1 = sigseedt(Nrandp(1:N1));
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal));
                    targetsig1 = dattarget(Nrandp(1:N1),:);
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:);
                    covt1 = covt(Nrandp(1:N1),:);
                    covt2 = covt(Nrandp(N1+1:Ntotal),:);
                    [r1 p1] = partialcorr(seedsig1,targetsig1,covt1);
                    [r2 p2] = partialcorr(seedsig2,targetsig2,covt2);
                    deltr(iperm,:) = r1-r2;
                    Rtemp1(iperm,:) = r1; 
                    Rtemp2(iperm,:) = r2;
                end
                save([Outdir,filesep,'ExtRval_',num2str(i),'.mat'],'Rtemp1','Rtemp2');
                [mu,sig,~,sigci] = normfit(deltr);
                P = normcdf(DELTR(i,:),mu,sig);
                Pval(i,:) = P;
            end
            INDEXS = computeval1.indexstarget;
            vmat = computeval1.vtarget;
            for i = 1:NROI
                if i<10
                    OutPname = fullfile(Outdir,['Pval_ROI0000',num2str(i),'.nii']);
                elseif i<100
                    OutPname = fullfile(Outdir,['Pval_ROI000',num2str(i),'.nii']);
                elseif i<1000
                    OutPname = fullfile(Outdir,['Pval_ROI00',num2str(i),'.nii']);
                else
                    error('too many ROIs');
                end
                Ptemp = zeros(vmat.dim);
                Ptemp(INDEXS) = Pval(i,:);
                DynamicBC_write_NIFTI(Ptemp,vmat,OutPname);
            end
        end
    end
else % interaction
    
end
% %% add in later version for the number of connect difference
% for i = 1:NROI
%     
% end

end