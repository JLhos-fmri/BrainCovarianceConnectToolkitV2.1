function AS_WTA_Stat_main(Parameter)
% Face a problem for the number of connects
% In later version, we would like to deal with it for the small ROM
% computer.
Outdir = Parameter.Outdir;
Indir1 = Parameter.Input1;
Indir2 = Parameter.Input2;
Incalmat1 = load(fullfile(Indir1,'SetUpparameter.mat')); 
Incalmat2 = load(fullfile(Indir2,'SetUpparameter.mat'));
COVcond1 = Incalmat1.Parameter.covs; % g1 Э�������
COVcond2 = Incalmat2.Parameter.covs; % g2 Э�������
Calmethod1 = Incalmat1.Parameter.methodused; % G1 ���㷽�� 1 pearson 2 partialcorr
Calmethod2 = Incalmat2.Parameter.methodused; % G2 ���㷽�� 1 pearson 2 partialcorr 
if COVcond1~=COVcond2
    error('The two groups used different cov conditions');
end
if Calmethod1~=Calmethod2
    error('The tow groups used different calculate methods');
end
COVCOND = COVcond1;
PermMethod = Parameter.permtype; % ͳ�Ʒ����� 1 permutation 2 interaction 

computeval1 = load(fullfile(Indir1,'computeval.mat'));
computeval2 = load(fullfile(Indir2,'computeval.mat'));
OrigR1 = computeval1.r; % G1 no scat��rֵ
OrigR2 = computeval2.r; % G2 no scat��rֵ
DELTR = OrigR1-OrigR2; % ����ֵ����Ҫ����permutation
SIG1 = load(fullfile(Indir1,'Signals.mat')); % ��ȡԭʼ�ź�
SIG2 = load(fullfile(Indir2,'Signals.mat')); % ��ȡԭʼ�ź�
datTarget1 = SIG1.datTarget; % G1 target�źţ�ԭʼ��
datTarget2 = SIG2.datTarget; % G2 target�źţ�ԭʼ��
datSeed1 = SIG1.datSeed; % G1 seed�źţ�ԭʼ��
datSeed2 = SIG2.datSeed; % G2 seed�źţ�ԭʼ��
NROI = size(datSeed1,1);
if PermMethod==1 % permutation 
    Permnum = Parameter.permnum; % permutation��Ŀ
    N1 = size(datTarget1,2); % G1 ��������ʱ��㳤�ȣ�
    N2 = size(datTarget2,2); % G2 ��������ʱ��㳤�ȣ�
    Ntotal = N1+N2; % ���屻����Ŀ��ʱ��㳤�ȣ�
    dattarget = [datTarget1,datTarget2]'; % ƴ�ӵ�target����
    if COVCOND==0 % no cov 
        if Calmethod1==1 % pearson & nocov
            for i = 1:NROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]'; % ��i��ROI��seed�ź�ƴ��
                deltr = zeros(Permnum,size(datTarget1,1)); % Ԥ��ֵdelta-pֵ����
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                parfor iperm = 1:Permnum % ����û�
                    Nrandp = randperm(Ntotal); % ����
                    seedsig1 = sigseedt(Nrandp(1:N1));  % ǰN1��Seed
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal)); % ��N2��Seed
                    targetsig1 = dattarget(Nrandp(1:N1),:);  % ǰN1��target
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:); % ��N2��target
                    [r1 p1] = corr(seedsig1,targetsig1); %��ط���
                    [r2 p2] = corr(seedsig2,targetsig2); %��ط���
                    deltr(iperm,:) = r1-r2; % delta-pֵ¼��
                    Rtemp1(iperm,:) = r1; 
                    Rtemp2(iperm,:) = r2;
                end
                save([Outdir,filesep,'ExtRval_',num2str(i),'.mat'],'Rtemp1','Rtemp2');
                [mu,sig,~,sigci] = normfit(deltr); %��̬�ֲ�
                P = normcdf(DELTR(i,:),mu,sig); % ����pֵ��<0.05ΪG1<G2��>0.95ΪG1>G2)
                Pval(i,:) = P; % ¼��pֵ
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
            for i = 1:NROI % ���ROI
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]';  % ƴ�����ӵ�
                deltr = zeros(Permnum,size(datTarget1,1));  % ����r�������
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                INDUSED = 1:NROI;INDUSED(i) = []; % ����Э����
                covt = [datSeed1(INDUSED,:),datSeed2(INDUSED,:)]'; % Э����ƴ��
                parfor iperm = 1:Permnum % ����û�
                    Nrandp = randperm(Ntotal); % �������
                    seedsig1 = sigseedt(Nrandp(1:N1)); % ǰN1��Seed
                    seedsig2 = sigseedt(Nrandp(N1+1:Ntotal)); % ��N2��seed
                    targetsig1 = dattarget(Nrandp(1:N1),:); % ǰN1��target
                    targetsig2 = dattarget(Nrandp(N1+1:Ntotal),:); % ��N2��target
                    covt1 = covt(Nrandp(1:N1),:); % ǰN1��COV
                    covt2 = covt(Nrandp(N1+1:Ntotal),:); % ��N2��COV
                    [r1 p1] = partialcorr(seedsig1,targetsig1,covt1); % ƫ���
                    [r2 p2] = partialcorr(seedsig2,targetsig2,covt2); % ƫ���
                    deltr(iperm,:) = r1-r2; % rֵ����¼��
                    Rtemp1(iperm,:) = r1; 
                    Rtemp2(iperm,:) = r2;
                end
                save([Outdir,filesep,'ExtRval_',num2str(i),'.mat'],'Rtemp1','Rtemp2');
                [mu,sig,~,sigci] = normfit(deltr); % ��̬�ֲ�
                P = normcdf(DELTR(i,:),mu,sig); % pֵ����
                Pval(i,:) = P; %pֵ¼��
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
                sigseedt = [datSeed1(i,:),datSeed2(i,:)]'; % ��i��ROI
                deltr = zeros(Permnum,size(datTarget1,1)); % ����delta-r
                Rtemp1 = zeros(Permnum,size(datTarget1,1));  
                Rtemp2 = zeros(Permnum,size(datTarget1,1));  
                covt = [SIG1.COV;SIG2.COV]; % Э����
                parfor iperm = 1:Permnum % ����û�
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
                covt = [covtCOV,[datSeed1(INDUSED,:),datSeed2(INDUSED,:)]']; % cov��COV������seed�źŹ���
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