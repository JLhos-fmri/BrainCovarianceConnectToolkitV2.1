function BCCT_SCN_stat_mainfunc(Parameter)
Outputdir = Parameter.Outputdir;
Inputdir1 = Parameter.Inputdir1;
Inputdir2 = Parameter.Inputdir2;

RealComp1 = load(fullfile(Inputdir1,'RealCompPara.mat'));
RealComp2 = load(fullfile(Inputdir2,'RealCompPara.mat'));
if ~strcmp(RealComp1.RealCompPara.mod,RealComp2.RealCompPara.mod)
    error('different SCN types in two input directory');
end

if strcmp(RealComp1.RealCompPara.mod,'volmap')
    load(fullfile(Inputdir1,'ROIsignal.mat'));
    ROIsig_g1 = ROIsignals;
    load(fullfile(Inputdir2,'ROIsignal.mat'));
    ROIsig_g2 = ROIsignals;
    load(fullfile(Inputdir1,'maskedSignal'));
    maskedsignals_g1 = maskedSignal;
    load(fullfile(Inputdir2,'maskedSignal'));
    maskedsignals_g2 = maskedSignal;
    vmask = RealComp1.RealCompPara.V(1);
    dims = RealComp1.RealCompPara.V(1).dim;
    DATMASK = RealComp1.RealCompPara.DATMASK;
    if Parameter.Intlab==1
        glab = 1;
        for i = 1:size(ROIsig_g1,1)
            Group{glab,1} = 'G1';
            glab = glab+1;
        end
        for i = 1:size(ROIsig_g2,1)
            Group{glab,1} = 'G2';
            glab = glab+1;
        end
        GROUP = term(Group);
        for i = 1:size(ROIsig_g1,2)
            if i<10
                Outfilenametemp = fullfile(Outputdir,['T_ROI00000',num2str(i),'.nii']);
                OutfilenametempZ = fullfile(Outputdir,['Z_ROI00000',num2str(i),'.nii']);
                OutfilenametempP = fullfile(Outputdir,['P_ROI00000',num2str(i),'.nii']);
            elseif i<100
                Outfilenametemp = fullfile(Outputdir,['T_ROI0000',num2str(i),'.nii']);
                OutfilenametempZ = fullfile(Outputdir,['Z_ROI0000',num2str(i),'.nii']);
                OutfilenametempP = fullfile(Outputdir,['P_ROI0000',num2str(i),'.nii']);
            elseif i<1000
                Outfilenametemp = fullfile(Outputdir,['T_ROI000',num2str(i),'.nii']);
                OutfilenametempZ = fullfile(Outputdir,['Z_ROI000',num2str(i),'.nii']);
                OutfilenametempP = fullfile(Outputdir,['P_ROI000',num2str(i),'.nii']);
            else
                Outfilenametemp = fullfile(Outputdir,['T_ROI00',num2str(i),'.nii']);
                OutfilenametempZ = fullfile(Outputdir,['Z_ROI00',num2str(i),'.nii']);
                OutfilenametempP = fullfile(Outputdir,['P_ROI00',num2str(i),'.nii']);
            end
            corrsig1 = ROIsig_g1(:,i);
            corrsig2 = ROIsig_g2(:,i);
            corrsig = [corrsig1;corrsig2];
            CORRSIG = term(corrsig);
            if RealComp1.RealCompPara.COVcond(1)
                COV1 = RealComp1.RealCompPara.COVs;
                COV2 = RealComp2.RealCompPara.COVs;
                covtotal = [COV1;COV2];
                COV = term(covtotal);
                mod = 1+GROUP+CORRSIG+GROUP*CORRSIG+COV;
                contrast = [0,0,0,0,1,-1,zeros(1,size(COV,2))];
            else
                mod = 1+GROUP+CORRSIG+GROUP*CORRSIG;
                contrast = [0,0,0,0,1,-1];
            end
            for j = 1:size(maskedSignal,2)
                try
                    Y = [maskedsignals_g1(:,j);maskedsignals_g2(:,j)];
                    slm = SurfStatLinMod( Y, mod );
                    slmt = SurfStatT(slm,contrast);
                    T(j) = slmt.t;
                    [Z(j) P(j)] = AS_TFRtoZ(slmt.t,'T',slmt.df,[]);
                catch
                    T(j) = 0;
                    Z(j) = 0;
                    P(j) = 0.5;
                end
            end
            Rmaps = zeros(dims);
            Rmaps(find(DATMASK)) = T;
            vmasknew = RealComp1.RealCompPara.V(1);
            vmasknew.descrip=sprintf('{T_[%.1f]}',slmt.df);
            DynamicBC_write_NIFTI(Rmaps,vmasknew,Outfilenametemp);
            Zmaps = zeros(dims);
            Zmaps(find(DATMASK)) = Z;
            DynamicBC_write_NIFTI(Zmaps,vmask,OutfilenametempZ);
            Pmaps = zeros(dims);
            Pmaps(find(DATMASK)) = P;
            DynamicBC_write_NIFTI(Pmaps,vmask,OutfilenametempP);
        end
    end
    if Parameter.Permlab==1
        Permnum = Parameter.PermNum;
        for i = 1:size(ROIsig_g1,2)
            if i<10
                Outfilenametemp = fullfile(Outputdir,['PermP_ROI00000',num2str(i),'.nii']);
            elseif i<100
                Outfilenametemp = fullfile(Outputdir,['PermP_ROI0000',num2str(i),'.nii']);
            elseif i<1000
                Outfilenametemp = fullfile(Outputdir,['PermP_ROI000',num2str(i),'.nii']);
            else
                Outfilenametemp = fullfile(Outputdir,['PermP_ROI00',num2str(i),'.nii']);
            end
            
            corrsig1 = ROIsig_g1(:,i);
            corrsig2 = ROIsig_g2(:,i);
            corrsig = [corrsig1;corrsig2];
            G1num = length(corrsig1);
            G2num = length(corrsig2);
            Totalnum = G1num+G2num;
            MASKEDSIGNALS = [maskedsignals_g1;maskedsignals_g2];
            
            if RealComp1.RealCompPara.COVcond(1)
                COV1 = RealComp1.RealCompPara.COVs;
                COV2 = RealComp2.RealCompPara.COVs;
                covtotal = [COV1;COV2];
                for iperm = 1:Permnum
                    randord = randperm(Totalnum);
                    g1sig = corrsig(randord(1:G1num),:);
                    g1maskedsig = MASKEDSIGNALS(randord(1:G1num),:);
                    g2sig = corrsig(randord(G1num+1:Totalnum),:);
                    g2maskedsig = MASKEDSIGNALS(randord(G1num+1:Totalnum),:);
                    cov1 = covtotal(randord(1:G1num),:);
                    cov2 = covtotal(randord(G1num+1:Totalnum),:);
                    [r1_perm(iperm,:) P] = partialcorr(g1maskedsig,g1sig,cov1);
                    [r2_perm(iperm,:) P] = partialcorr(g2maskedsig,g2sig,cov2);
                end
                [R1,P1] = partialcorr(maskedsignals_g1,corrsig1,COV1);
                [R2,P2] = partialcorr(maskedsignals_g2,corrsig2,COV2);
                [mu,sig,~,sigci] = normfit(r1_perm-r2_perm);
                P_map1 = normcdf((R1-R2)',mu,sig);
            else
                %                [R P] = corr(maskedSignal,ROIsignals(:,i));
                %                DF_E = size(ROIsignals,1)-2;
                %                [Z P2] = AS_TFRtoZ(R,'R',DF_E,[]);
                for iperm = 1:Permnum
                    randord = randperm(Totalnum);
                    g1sig = corrsig(randord(1:G1num),:);
                    g1maskedsig = MASKEDSIGNALS(randord(1:G1num),:);
                    g2sig = corrsig(randord(G1num+1:Totalnum),:);
                    g2maskedsig = MASKEDSIGNALS(randord(G1num+1:Totalnum),:);
%                     cov1 = covtotal(randord(1:G1num),:);
%                     cov2 = covtotal(randord(G1num+1:Totalnum),:);
                    [r1_perm(iperm,:) P] = corr(g1maskedsig,g1sig);
                    [r2_perm(iperm,:) P] = corr(g2maskedsig,g2sig);
                end
                [R1,P1] = corr(maskedsignals_g1,corrsig1);
                [R2,P2] = corr(maskedsignals_g2,corrsig2);
                [mu,sig,~,sigci] = normfit(r1_perm-r2_perm);
                P_map1 = normcdf((R1-R2)',mu,sig);
            end
            Pmaps = zeros(dims);
            Pmaps(find(DATMASK)) = P_map1;
            DynamicBC_write_NIFTI(Pmaps,vmask,Outfilenametemp);
        end
    end
    
elseif strcmp(RealComp1.RealCompPara.mod,'volroi')    
    load(fullfile(Inputdir1,'ROIsignal.mat'));
    ROIsig_g1 = ROIsignals;
    load(fullfile(Inputdir2,'ROIsignal.mat'));
    ROIsig_g2 = ROIsignals;
    
    
    if Parameter.Intlab==1        
        glab = 1;
        for i = 1:size(ROIsig_g1,1)
            Group{glab,1} = 'G1';
            glab = glab+1;
        end
        for i = 1:size(ROIsig_g2,1)
            Group{glab,1} = 'G2';
            glab = glab+1;
        end
        GROUP = term(Group);
        if RealComp1.RealCompPara.Partinfo==1 % partial correlation
            for i = 1:size(ROIsig_g1,2)
                corrsig1 = ROIsig_g1(:,i);
                corrsig2 = ROIsig_g2(:,i);
                corrsig = [corrsig1;corrsig2];
                CORRSIG = term(corrsig);
                if RealComp1.RealCompPara.COVcond(1)      % with cov              
                    COV1 = RealComp1.RealCompPara.COVs;
                    COV2 = RealComp2.RealCompPara.COVs;
                    covtotal = [COV1;COV2];
%                     COV = term(covtotal);
                    for j = 1:size(ROIsig_g1,2)
                        if i==j
                            T(i,j) = 0;
                            P(i,j) = 1;
                            Z(i,j) = 0;
                        else
                            leftnum = 1:size(ROIsig_g1,2);
                            leftnum([i,j]) = [];
                            covtotal1 = [covtotal,[ROIsig_g1(:,leftnum);ROIsig_g2(:,leftnum)]];
                            COV = term(covtotal1);
                            mod = 1+GROUP+CORRSIG+GROUP*CORRSIG+COV;
                            contrast = [0,0,0,0,1,-1,zeros(1,size(COV,2))];
                            try
                                Y = [ROIsig_g1(:,j);ROIsig_g2(:,j)];
                                slm = SurfStatLinMod( Y, mod );
                                slmt = SurfStatT(slm,contrast);
                                T(i,j) = slmt.t;
                                [Z(i,j) P(i,j)] = AS_TFRtoZ(slmt.t,'T',slmt.df,[]);
                            catch
                                T(i,j) = 0;
                                Z(i,j) = 0;
                                P(i,j) = 0.5;
                            end
                        end
                    end                    
                else
                    for j = 1:size(ROIsig_g1,2)
                        if i==j
                            T(i,j) = 0;
                            P(i,j) = 1;
                            Z(i,j) = 0;
                        else
                            leftnum = 1:size(ROIsig_g1,2)
                            leftnum([i,j]) = [];
                            covtotal1 = [ROIsig_g1(:,leftnum);ROIsig_g2(:,left)];
                            COV = term(covtotal1);
                            mod = 1+GROUP+CORRSIG+GROUP*CORRSIG+COV;
                            contrast = [0,0,0,0,1,-1,zeros(1,size(COV,2))];
                            try
                                Y = [ROIsig_g1(:,j);ROIsig_g2(:,j)];
                                slm = SurfStatLinMod( Y, mod );
                                slmt = SurfStatT(slm,contrast);
                                T(i,j) = slmt.t;
                                [Z(i,j) P(i,j)] = AS_TFRtoZ(slmt.t,'T',slmt.df,[]);
                            catch
                                T(i,j) = 0;
                                Z(i,j) = 0;
                                P(i,j) = 0.5;
                            end
                        end
                    end  
                end
            end
        else
            for i = 1:size(ROIsig_g1,2)
                corrsig1 = ROIsig_g1(:,i);
                corrsig2 = ROIsig_g2(:,i);
                corrsig = [corrsig1;corrsig2];
                CORRSIG = term(corrsig);
                if RealComp1.RealCompPara.COVcond(1)      % with cov              
                    COV1 = RealComp1.RealCompPara.COVs;
                    COV2 = RealComp2.RealCompPara.COVs;
                    covtotal = [COV1;COV2];
                    COV = term(covtotal);
                    mod = 1+GROUP+CORRSIG+GROUP*CORRSIG+COV;
                    contrast = [0,0,0,0,1,-1,zeros(1,size(COV,2))];
                else
                    mod = 1+GROUP+CORRSIG+GROUP*CORRSIG;
                    contrast = [0,0,0,0,1,-1];
                end
                try
                    Y = [ROIsig_g1;ROIsig_g2];
                    slm = SurfStatLinMod( Y, mod );
                    slmt = SurfStatT(slm,contrast);
                    T(i,:) = slmt.t;
                    [Z(i,:) P(i,:)] = AS_TFRtoZ(slmt.t,'T',slmt.df,[]);
                catch
                    T(i,:) = 0;
                    Z(i,:) = 0;
                    P(i,:) = 0.5;
                end
            end
        end
        save(fullfile(Outputdir,'InterCompRes.mat'),'T','P','Z')
    end
    if Parameter.Permlab==1
        Permnum = Parameter.PermNum;
        G1num = size(ROIsig_g1,1);
        G2num = size(ROIsig_g2,1);
        Totalnum = G1num+G2num;
        totalsig = [ROIsig_g1;ROIsig_g2];
        if RealComp1.RealCompPara.Partinfo==1 % partial correlation
            if RealComp1.RealCompPara.COVcond(1)      % with cov        
                COV1 = RealComp1.RealCompPara.COVs;
                COV2 = RealComp2.RealCompPara.COVs;
                covtotal = [COV1;COV2];
                
                for iperm = 1:Permnum
                    randord = randperm(Totalnum);
                    G1ord = randord(1:G1num);
                    G2ord = randord(G1num+1:Totalnum);
                    sig_lab = 1;
                    for i = 1:size(totalsig,2)-1
                        for j = i+1:size(totalsig,2)                            
                            leftnum = 1:size(ROIsig_g1,2);
                            leftnum([i,j]) = [];
                            COVUsed = [covtotal,totalsig(:,leftnum)];
                            [r1 p1] = partialcorr(totalsig(G1ord,i),totalsig(G1ord,j),COVUsed(G1ord,:));
                            [r2 p2] = partialcorr(totalsig(G2ord,i),totalsig(G2ord,j),COVUsed(G2ord,:));
                            rdiff(iperm,sig_lab) = r1-r2;
                            sig_lab = sig_lab+1;
                        end
                    end
                end
                R_g1 = load(fullfile(Inputdir1,'R_Pres.mat'));
                R_g2 = load(fullfile(Inputdir2,'R_Pres.mat'));
                sig_lab = 1;
                P_mat = ones(size(totalsig,2))*0.5;
                for i = 1:size(totalsig,2)-1
                    for j = i+1:size(totalsig,2)
                        Rdiff = R_g1.R(i,j)-R_g2.R(i,j);
                        [mu,sig,~,sigci] = normfit(rdiff(:,sig_lab));
                        P_mat(i,j) = normcdf(Rdiff,mu,sig);
                        P_mat(j,i) = P_mat(i,j);
                        sig_lab = sig_lab+1;
                    end
                end
                
            else
                for iperm = 1:Permnum
                    randord = randperm(Totalnum);
                    G1ord = randord(1:G1num);
                    G2ord = randord(G1num+1:Totalnum);
                    sig_lab = 1;
                    for i = 1:size(totalsig,2)-1
                        for j = i+1:size(totalsig,2)                            
                            leftnum = 1:size(ROIsig_g1,2);
                            leftnum([i,j]) = [];
                            COVUsed = [totalsig(:,leftnum)];
                            [r1 p1] = partialcorr(totalsig(G1ord,i),totalsig(G1ord,j),COVUsed(G1ord,:));
                            [r2 p2] = partialcorr(totalsig(G2ord,i),totalsig(G2ord,j),COVUsed(G1ord,:));
                            rdiff(iperm,sig_lab) = r1-r2;
                            sig_lab = sig_lab+1;
                        end
                    end
                end
                
                R_g1 = load(fullfile(Inputdir1,'R_Pres.mat'));
                R_g2 = load(fullfile(Inputdir2,'R_Pres.mat'));
                sig_lab = 1;
                P_mat = ones(size(totalsig,2))*0.5;
                for i = 1:size(totalsig,2)-1
                    for j = i+1:size(totalsig,2)
                        Rdiff = R_g1.R(i,j)-R_g2.R(i,j);
                        [mu,sig,~,sigci] = normfit(rdiff(:,sig_lab));
                        P_mat(i,j) = normcdf(Rdiff,mu,sig);
                        P_mat(j,i) = P_mat(i,j);
                        sig_lab = sig_lab+1;
                    end
                end
            end
        else
            if RealComp1.RealCompPara.COVcond(1)      % with cov    
                COV1 = RealComp1.RealCompPara.COVs;
                COV2 = RealComp2.RealCompPara.COVs;
                covtotal = [COV1;COV2];
                for iperm = 1:Permnum
                    randord = randperm(Totalnum);
                    G1ord = randord(1:G1num);
                    G2ord = randord(G1num+1:Totalnum);
                    [r1 p] = partialcorr(totalsig(G1ord,:),totalsig(G1ord,:),covtotal(G1ord,:));
                    [r2 p] = partialcorr(totalsig(G2ord,:),totalsig(G2ord,:),covtotal(G2ord,:));
                    difr = r1-r2;
                    rdiff(iperm,:) = reshape(difr,1,size(difr,1)*size(difr,2));
                end
                
                R_g1 = load(fullfile(Inputdir1,'R_Pres.mat'));
                R_g2 = load(fullfile(Inputdir2,'R_Pres.mat'));
%                 P_mat = ones(size(totalsig,2))*0.5;
                Rdiff = reshape(R_g1.R-R_g2.R,1,size(R_g1.R,1)*size(R_g1.R,2));
                [mu,sig] = normfit(rdiff);
                P_mat1 = normcdf(Rdiff,mu,sig);
                P_mat = reshape(P_mat1,size(R_g1.R,1),size(R_g1.R,2));
            else
                for iperm = 1:Permnum
                    randord = randperm(Totalnum);
                    G1ord = randord(1:G1num);
                    G2ord = randord(G1num+1:Totalnum);
                    [r1 p] = corr(totalsig(G1ord,:),totalsig(G1ord,:));
                    [r2 p] = corr(totalsig(G2ord,:),totalsig(G2ord,:));
                    difr = r1-r2;
                    rdiff(iperm,:) = reshape(difr,1,size(difr,1)*size(difr,2));
                end
                R_g1 = load(fullfile(Inputdir1,'R_Pres.mat'));
                R_g2 = load(fullfile(Inputdir2,'R_Pres.mat'));
%                 P_mat = ones(size(totalsig,2))*0.5;
                Rdiff = reshape(R_g1.R-R_g2.R,1,size(R_g1.R,1)*size(R_g1.R,2));
                [mu,sig] = normfit(rdiff);
                P_mat1 = normcdf(Rdiff,mu,sig);
                P_mat = reshape(P_mat1,size(R_g1.R,1),size(R_g1.R,2));
            end
        end
        
        save(fullfile(Outputdir,'PermCompRes.mat'),'P_mat')
    end
elseif strcmp(RealComp1.RealCompPara.mod,'surfmap')
    
    load(fullfile(Inputdir1,'ROIsignal.mat'));
    ROIsig_g1 = ROIsignals;
    load(fullfile(Inputdir2,'ROIsignal.mat'));
    ROIsig_g2 = ROIsignals;
    load(fullfile(Inputdir1,'maskedSignal'));
    maskedsignals_g1 = maskedSignal;
    load(fullfile(Inputdir2,'maskedSignal'));
    maskedsignals_g2 = maskedSignal;
    MASK = RealComp1.RealCompPara.MASK;
    Ml = RealComp1.RealCompPara.Ml;
    Mr = RealComp1.RealCompPara.Mr;
    N1size = RealComp1.RealCompPara.N1size;
    N2size = RealComp1.RealCompPara.N2size;
%     RealCompPara.Ml = Ml;
%     RealCompPara.Mr = Mr;
%     RealCompPara.N1size = N1size;
%     RealCompPara.N2size = N2size;
    if Parameter.Intlab==1
        glab = 1;
        for i = 1:size(ROIsig_g1,1)
            Group{glab,1} = 'G1';
            glab = glab+1;
        end
        for i = 1:size(ROIsig_g2,1)
            Group{glab,1} = 'G2';
            glab = glab+1;
        end
        GROUP = term(Group);
%         save temp
        for i = 1:size(ROIsig_g1,2)
            if i<10
                Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI00000',num2str(i),'.mgh']);
                OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI00000',num2str(i),'.mgh']);
                OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI00000',num2str(i),'.mgh']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI00000',num2str(i),'.mgh']);
                OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI00000',num2str(i),'.mgh']);
                OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI00000',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['mat_ROI00000',num2str(i),'.mat']);
                
            elseif i<100
                Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI0000',num2str(i),'.mgh']);
                OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI0000',num2str(i),'.mgh']);
                OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI0000',num2str(i),'.mgh']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI0000',num2str(i),'.mgh']);
                OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI0000',num2str(i),'.mgh']);
                OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI0000',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['mat_ROI0000',num2str(i),'.mat']);
            elseif i<1000
                Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI000',num2str(i),'.mgh']);
                OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI000',num2str(i),'.mgh']);
                OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI000',num2str(i),'.mgh']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI000',num2str(i),'.mgh']);
                OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI000',num2str(i),'.mgh']);
                OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI000',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['mat_ROI000',num2str(i),'.mat']);
            else
                Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI00',num2str(i),'.mgh']);
                OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI00',num2str(i),'.mgh']);
                OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI00',num2str(i),'.mgh']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI00',num2str(i),'.mgh']);
                OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI00',num2str(i),'.mgh']);
                OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI00',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['mat_ROI00',num2str(i),'.mat']);
            end
            corrsig1 = ROIsig_g1(:,i);
            corrsig2 = ROIsig_g2(:,i);
            corrsig = [corrsig1;corrsig2];
            CORRSIG = term(corrsig);
            if RealComp1.RealCompPara.COVcond(1)
                COV1 = RealComp1.RealCompPara.COVs;
                COV2 = RealComp2.RealCompPara.COVs;
                covtotal = [COV1;COV2];
                COV = term(covtotal);
                mod = 1+GROUP+CORRSIG+GROUP*CORRSIG+COV;
                contrast = [0,0,0,0,1,-1,zeros(1,size(COV,2))];
            else
                mod = 1+GROUP+CORRSIG+GROUP*CORRSIG;
                contrast = [0,0,0,0,1,-1];
            end
            for j = 1:size(maskedSignal,2)
                try
                    Y = [maskedsignals_g1(:,j);maskedsignals_g2(:,j)];
                    slm = SurfStatLinMod( Y, mod );
                    slmt = SurfStatT(slm,contrast);
                    T(j) = slmt.t;
                    [Z(j) P(j)] = AS_TFRtoZ(slmt.t,'T',slmt.df,[]);
                catch
                    T(j) = 0;
                    Z(j) = 0;
                    P(j) = 0.5;
                end
            end
            Tmaps = zeros(size(MASK));
            Tmaps(find(MASK)) = T;
            
            save_mgh(Tmaps(1:N1size),Outfilenametemp_L,Ml);
            save_mgh(Tmaps(1+N1size:N1size+N2size),Outfilenametemp_R,Mr);
%             SurfStatWriteData(Outfilenametemp,Tmaps);
            Zmaps = zeros(size(MASK));
            Zmaps(find(MASK)) = Z;
            
            save_mgh(Zmaps(1:N1size),OutfilenametempZ_L,Ml);
            save_mgh(Zmaps(1+N1size:N1size+N2size),OutfilenametempZ_R,Mr);
%             SurfStatWriteData(OutfilenametempZ,Zmaps);
            Pmaps = zeros(size(MASK));
            Pmaps(find(MASK)) = P;
            
            save_mgh(Pmaps(1:N1size),OutfilenametempP_L,Ml);
            save_mgh(Pmaps(1+N1size:N1size+N2size),OutfilenametempP_R,Mr);
%             SurfStatWriteData(OutfilenametempP,Pmaps);
            save(OutfilenametempMat,'Tmaps','Zmaps','Pmaps')
        end
    end
    if Parameter.Permlab==1
        Permnum = Parameter.PermNum;
        for i = 1:size(ROIsig_g1,2)
            if i<10
                Outfilenametemp_L = fullfile(Outputdir,['lh.PermP_ROI00000',num2str(i),'.mgh']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.PermP_ROI00000',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['PermP_ROI00000',num2str(i),'.mat']);
            elseif i<100
                Outfilenametemp_L = fullfile(Outputdir,['lh.PermP_ROI0000',num2str(i),'.mgh']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.PermP_ROI0000',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['PermP_ROI0000',num2str(i),'.mat']);
            elseif i<1000
                Outfilenametemp_L = fullfile(Outputdir,['lh.PermP_ROI000',num2str(i),'.mgh']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.PermP_ROI000',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['PermP_ROI000',num2str(i),'.mat']);
            else
                Outfilenametemp_L = fullfile(Outputdir,['lh.PermP_ROI00',num2str(i),'.mgh']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.PermP_ROI00',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['PermP_ROI00',num2str(i),'.mat']);
            end
            
            corrsig1 = ROIsig_g1(:,i);
            corrsig2 = ROIsig_g2(:,i);
            corrsig = [corrsig1;corrsig2];
            G1num = length(corrsig1);
            G2num = length(corrsig2);
            Totalnum = G1num+G2num;
            MASKEDSIGNALS = [maskedsignals_g1;maskedsignals_g2];
            
            if RealComp1.RealCompPara.COVcond(1)
                COV1 = RealComp1.RealCompPara.COVs;
                COV2 = RealComp2.RealCompPara.COVs;
                covtotal = [COV1;COV2];
                for iperm = 1:Permnum
                    randord = randperm(Totalnum);
                    g1sig = corrsig(randord(1:G1num),:);
                    g1maskedsig = MASKEDSIGNALS(randord(1:G1num),:);
                    g2sig = corrsig(randord(G1num+1:Totalnum),:);
                    g2maskedsig = MASKEDSIGNALS(randord(G1num+1:Totalnum),:);
                    cov1 = covtotal(randord(1:G1num),:);
                    cov2 = covtotal(randord(G1num+1:Totalnum),:);
                    [r1_perm(iperm,:) P] = partialcorr(g1maskedsig,g1sig,cov1);
                    [r2_perm(iperm,:) P] = partialcorr(g2maskedsig,g2sig,cov2);
                end
                [R1,P1] = partialcorr(maskedsignals_g1,corrsig1,COV1);
                [R2,P2] = partialcorr(maskedsignals_g2,corrsig2,COV2);
                [mu,sig,~,sigci] = normfit(r1_perm-r2_perm);
                P_map1 = normcdf((R1-R2)',mu,sig);
            else
                %                [R P] = corr(maskedSignal,ROIsignals(:,i));
                %                DF_E = size(ROIsignals,1)-2;
                %                [Z P2] = AS_TFRtoZ(R,'R',DF_E,[]);
                for iperm = 1:Permnum
                    randord = randperm(Totalnum);
                    g1sig = corrsig(randord(1:G1num),:);
                    g1maskedsig = MASKEDSIGNALS(randord(1:G1num),:);
                    g2sig = corrsig(randord(G1num+1:Totalnum),:);
                    g2maskedsig = MASKEDSIGNALS(randord(G1num+1:Totalnum),:);
%                     cov1 = covtotal(randord(1:G1num),:);
%                     cov2 = covtotal(randord(G1num+1:Totalnum),:);
                    [r1_perm(iperm,:) P] = corr(g1maskedsig,g1sig);
                    [r2_perm(iperm,:) P] = corr(g2maskedsig,g2sig);
                end
                [R1,P1] = corr(maskedsignals_g1,corrsig1);
                [R2,P2] = corr(maskedsignals_g2,corrsig2);
                [mu,sig,~,sigci] = normfit(r1_perm-r2_perm);
                P_map1 = normcdf((R1-R2)',mu,sig);
            end
            Pmaps = zeros(size(MASK));
            Pmaps(find(MASK)) = P;
            
            save_mgh(Pmaps(1:N1size),Outfilenametemp_L,Ml);
            save_mgh(Pmaps(1+N1size:N1size+N2size),Outfilenametemp_R,Mr);
%             SurfStatWriteData(Outfilenametemp,Pmaps);
            save(OutfilenametempMat,'Pmaps')
        end
    end
elseif strcmp(RealComp1.RealCompPara.mod,'surfROI')
    disp('comming soon')
end
end