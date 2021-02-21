function BCCT_Modulation_stat_mainfunc(Parameter)

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
        modfactor1 = RealComp1.RealCompPara.Modufactor;
        modfactor2 = RealComp2.RealCompPara.Modufactor;
        modfactorall = [modfactor1;modfactor2];
        MODFACTOR = term(modfactorall);
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
                mod = 1+GROUP+CORRSIG+MODFACTOR+GROUP*CORRSIG+GROUP*MODFACTOR+CORRSIG*MODFACTOR+GROUP*MODFACTOR*CORRSIG+COV;
                contrast = [0,0,0,0,0,0,0,0,0,0,1,-1,zeros(1,size(COV,2))];
            else
                mod = 1+GROUP+CORRSIG+MODFACTOR+GROUP*CORRSIG+GROUP*MODFACTOR+CORRSIG*MODFACTOR+GROUP*MODFACTOR*CORRSIG;
                contrast = [0,0,0,0,0,0,0,0,0,0,1,-1];
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
    
elseif strcmp(RealComp1.RealCompPara.mod,'volroi')    
    load(fullfile(Inputdir1,'ROIsignal.mat'));
    ROIsig_g1 = ROIsignals;
    load(fullfile(Inputdir2,'ROIsignal.mat'));
    ROIsig_g2 = ROIsignals;
    
    modfactor1 = RealComp1.RealCompPara.Modufactor;
    modfactor2 = RealComp2.RealCompPara.Modufactor;
    modfactorall = [modfactor1;modfactor2];
    MODFACTOR = term(modfactorall);
    
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
%             for i = 1:size(ROIsig_g1,2)
%                 corrsig1 = ROIsig_g1(:,i);
%                 corrsig2 = ROIsig_g2(:,i);
%                 corrsig = [corrsig1;corrsig2];
%                 CORRSIG = term(corrsig);
%                 if RealComp1.RealCompPara.COVcond(1)      % with cov              
%                     COV1 = RealComp1.RealCompPara.COVs;
%                     COV2 = RealComp2.RealCompPara.COVs;
%                     covtotal = [COV1;COV2];
% %                     COV = term(covtotal);
%                     for j = 1:size(ROIsig_g1,2)
%                         if i==j
%                             T(i,j) = 0;
%                             P(i,j) = 1;
%                             Z(i,j) = 0;
%                         else
%                             leftnum = 1:size(ROIsig_g1,2);
%                             leftnum([i,j]) = [];
%                             covtotal1 = [covtotal,[ROIsig_g1(:,leftnum);ROIsig_g2(:,left)]];
%                             COV = term(covtotal1);
%                             mod = 1+GROUP+CORRSIG+GROUP*CORRSIG+COV;
%                             contrast = [0,0,0,0,1,-1,zeros(1,size(COV,2))];
%                             try
%                                 Y = [ROIsig_g1(:,j);ROIsig_g2(:,j)];
%                                 slm = SurfStatLinMod( Y, mod );
%                                 slmt = SurfStatT(slm,contrast);
%                                 T(i,j) = slmt.t;
%                                 [Z(i,j) P(i,j)] = AS_TFRtoZ(slmt.t,'T',slmt.df,[]);
%                             catch
%                                 T(i,j) = 0;
%                                 Z(i,j) = 0;
%                                 P(i,j) = 0.5;
%                             end
%                         end
%                     end                    
%                 else
%                     for j = 1:size(ROIsig_g1,2)
%                         if i==j
%                             T(i,j) = 0;
%                             P(i,j) = 1;
%                             Z(i,j) = 0;
%                         else
%                             leftnum = 1:size(ROIsig_g1,2)
%                             leftnum([i,j]) = [];
%                             covtotal1 = [ROIsig_g1(:,leftnum);ROIsig_g2(:,left)];
%                             COV = term(covtotal1);
%                             mod = 1+GROUP+CORRSIG+GROUP*CORRSIG+COV;
%                             contrast = [0,0,0,0,1,-1,zeros(1,size(COV,2))];
%                             try
%                                 Y = [ROIsig_g1(:,j);ROIsig_g2(:,j)];
%                                 slm = SurfStatLinMod( Y, mod );
%                                 slmt = SurfStatT(slm,contrast);
%                                 T(i,j) = slmt.t;
%                                 [Z(i,j) P(i,j)] = AS_TFRtoZ(slmt.t,'T',slmt.df,[]);
%                             catch
%                                 T(i,j) = 0;
%                                 Z(i,j) = 0;
%                                 P(i,j) = 0.5;
%                             end
%                         end
%                     end  
%                 end
%             end
        else
            for i = 1:size(ROIsig_g1,2)
                corrsig1 = ROIsig_g1(:,i);
                corrsig2 = ROIsig_g2(:,i);
                corrsig = [corrsig1;corrsig2];
                CORRSIG = term(corrsig);
                
                if RealComp1.RealCompPara.COVcond(1)  % with cov              
                    COV1 = RealComp1.RealCompPara.COVs;
                    COV2 = RealComp2.RealCompPara.COVs;
                    covtotal = [COV1;COV2];
                    COV = term(covtotal);
                    mod = 1+GROUP+CORRSIG+MODFACTOR+GROUP*CORRSIG+GROUP*MODFACTOR+CORRSIG*MODFACTOR+GROUP*MODFACTOR*CORRSIG+COV;
                    contrast = [0,0,0,0,0,0,0,0,0,0,1,-1,zeros(1,size(COV,2))];
                else
                    mod = 1+GROUP+CORRSIG+MODFACTOR+GROUP*CORRSIG+GROUP*MODFACTOR+CORRSIG*MODFACTOR+GROUP*MODFACTOR*CORRSIG;
                    contrast = [0,0,0,0,0,0,0,0,0,0,1,-1];
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
        save(fullfile(Outputdir,'ResultModfact.mat'),'T','Z','P');
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
        modfactor1 = RealComp1.RealCompPara.Modufactor;
        modfactor2 = RealComp2.RealCompPara.Modufactor;
        modfactorall = [modfactor1;modfactor2];
        MODFACTOR = term(modfactorall);
%         save temp
        for i = 1:size(ROIsig_g1,2)
            if i<10
                Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI00000',num2str(i),'.mgh']);
                OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI00000',num2str(i),'.mgh']);
                OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI00000',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['mat_ROI00000',num2str(i),'.mat']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI00000',num2str(i),'.mgh']);
                OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI00000',num2str(i),'.mgh']);
                OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI00000',num2str(i),'.mgh']);
                
            elseif i<100
                Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI0000',num2str(i),'.mgh']);
                OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI0000',num2str(i),'.mgh']);
                OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI0000',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['mat_ROI0000',num2str(i),'.mat']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI0000',num2str(i),'.mgh']);
                OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI0000',num2str(i),'.mgh']);
                OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI0000',num2str(i),'.mgh']);
            elseif i<1000
                Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI000',num2str(i),'.mgh']);
                OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI000',num2str(i),'.mgh']);
                OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI000',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['mat_ROI000',num2str(i),'.mat']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI000',num2str(i),'.mgh']);
                OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI000',num2str(i),'.mgh']);
                OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI000',num2str(i),'.mgh']);
            else
                Outfilenametemp_L = fullfile(Outputdir,['lh.T_ROI00',num2str(i),'.mgh']);
                OutfilenametempZ_L = fullfile(Outputdir,['lh.Z_ROI00',num2str(i),'.mgh']);
                OutfilenametempP_L = fullfile(Outputdir,['lh.P_ROI00',num2str(i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['mat_ROI00',num2str(i),'.mat']);
                Outfilenametemp_R = fullfile(Outputdir,['rh.T_ROI00',num2str(i),'.mgh']);
                OutfilenametempZ_R = fullfile(Outputdir,['rh.Z_ROI00',num2str(i),'.mgh']);
                OutfilenametempP_R = fullfile(Outputdir,['rh.P_ROI00',num2str(i),'.mgh']);
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
                mod = 1+GROUP+CORRSIG+MODFACTOR+GROUP*CORRSIG+GROUP*MODFACTOR+CORRSIG*MODFACTOR+GROUP*MODFACTOR*CORRSIG+COV;
                contrast = [0,0,0,0,0,0,0,0,0,0,1,-1,zeros(1,size(COV,2))];
            else
                mod = 1+GROUP+CORRSIG+MODFACTOR+GROUP*CORRSIG+GROUP*MODFACTOR+CORRSIG*MODFACTOR+GROUP*MODFACTOR*CORRSIG;
                contrast = [0,0,0,0,0,0,0,0,0,0,1,-1];
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
elseif strcmp(RealComp1.RealCompPara.mod,'surfROI')
    disp('comming soon')
end
end