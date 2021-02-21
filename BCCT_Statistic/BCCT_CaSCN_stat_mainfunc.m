function BCCT_CaSCN_stat_mainfunc(Parameter)
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
    MASK_G1 = DATMASK;
    maskedsignals_g1 = maskedSignal;
    load(fullfile(Inputdir2,'maskedSignal'));
    MASK_G2 = DATMASK;
    maskedsignals_g2 = maskedSignal;
    vmask = RealComp1.RealCompPara.V(1);
    dims = RealComp1.RealCompPara.V(1).dim;
    DATMASK = MASK_G1;
%     DATMASK = RealComp1.RealCompPara.DATMASK;
    if Parameter.Permlab==1
        Permnum = Parameter.PermNum;
        for i = 1:size(ROIsig_g1,2)    
            Res1orig = BCCT_Perm_CaSCNMap_Cal(RealComp1.RealCompPara,maskedsignals_g1,ROIsig_g1(:,i));
            Res2orig = BCCT_Perm_CaSCNMap_Cal(RealComp2.RealCompPara,maskedsignals_g2,ROIsig_g2(:,i));
            
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
            end
            for iperm = 1:Permnum
                randord = randperm(Totalnum);
                Ord1 = randord(1:G1num);
                Ord2 = randord(G1num+1:Totalnum);
                corrsig1perm = corrsig(Ord1,:);
                corrsig2perm = corrsig(Ord2,:);
                maskedsignal1perm = MASKEDSIGNALS(Ord1,:);
                maskedsignal2perm = MASKEDSIGNALS(Ord2,:);
                RealCompPara1temp = RealComp1.RealCompPara;
                RealCompPara2temp = RealComp2.RealCompPara;
                if RealComp1.RealCompPara.COVcond(1)
                    RealCompPara1temp.COVs = covtotal(Ord1,:);
                    RealCompPara2temp.COVs = covtotal(Ord2,:);
                end
                RealCompPara1temp.PermMark=0;
                RealCompPara2temp.PermMark=0;
                
                Res1(iperm) = BCCT_Perm_CaSCNMap_Cal(RealCompPara1temp,maskedsignal1perm,corrsig1perm(:,i));
                Res2(iperm) = BCCT_Perm_CaSCNMap_Cal(RealCompPara2temp,maskedsignal2perm,corrsig2perm(:,i));

            end
            if RealComp1.RealCompPara.Calmethod1 % res
                
                Outfilenametemp1 = fullfile(Outputdir,['PermP_x2y_ROI',sprintf('%05d',i),'.nii']);
                Outfilenametemp2 = fullfile(Outputdir,['PermP_y2x_ROI',sprintf('%05d',i),'.nii']);
                Outfilenametemp3 = fullfile(Outputdir,['PermP_x2ytrans_ROI',sprintf('%05d',i),'.nii']);
                Outfilenametemp4 = fullfile(Outputdir,['PermP_y2xtrans_ROI',sprintf('%05d',i),'.nii']);
                Outfilenametemp5 = fullfile(Outputdir,['PermP_Netx2y_ROI',sprintf('%05d',i),'.nii']);
                
                
                for iperm = 1:Permnum
                    G1_Map1temp(iperm,:) = Res1(iperm).res.ResultMap1;
                    G1_Map2temp(iperm,:) = Res1(iperm).res.ResultMap2;
                    G1_Map3temp(iperm,:) = Res1(iperm).res.ResultMap3;
                    G1_Map4temp(iperm,:) = Res1(iperm).res.ResultMap4;
                    G1_Map5temp(iperm,:) = Res1(iperm).res.ResultMap5;
                    G2_Map1temp(iperm,:) = Res2(iperm).res.ResultMap1;
                    G2_Map2temp(iperm,:) = Res2(iperm).res.ResultMap2;
                    G2_Map3temp(iperm,:) = Res2(iperm).res.ResultMap3;
                    G2_Map4temp(iperm,:) = Res2(iperm).res.ResultMap4;
                    G2_Map5temp(iperm,:) = Res2(iperm).res.ResultMap5;
                end
                [mu sig] = normfit(G1_Map1temp-G2_Map1temp);
                P_map1 = normcdf(Res1orig.res.ResultMap1-Res2orig.res.ResultMap1,mu,sig);
                DAT = zeros(dims);
                DAT(find(DATMASK)) = P_map1;
                DynamicBC_write_NIFTI(DAT,vmask,Outfilenametemp1)
                
                [mu sig] = normfit(G1_Map2temp-G2_Map2temp);
                P_map2 = normcdf(Res1orig.res.ResultMap2-Res2orig.res.ResultMap2,mu,sig);
                DAT = zeros(dims);
                DAT(find(DATMASK)) = P_map2;
                DynamicBC_write_NIFTI(DAT,vmask,Outfilenametemp2)
                
                [mu sig] = normfit(G1_Map3temp-G2_Map3temp);
                P_map3 = normcdf(Res1orig.res.ResultMap3-Res2orig.res.ResultMap3,mu,sig);
                DAT = zeros(dims);
                DAT(find(DATMASK)) = P_map3;
                DynamicBC_write_NIFTI(DAT,vmask,Outfilenametemp3)
                
                [mu sig] = normfit(G1_Map4temp-G2_Map4temp);
                P_map4 = normcdf(Res1orig.res.ResultMap4-Res2orig.res.ResultMap4,mu,sig);
                DAT = zeros(dims);
                DAT(find(DATMASK)) = P_map4;
                DynamicBC_write_NIFTI(DAT,vmask,Outfilenametemp4)
                
                [mu sig] = normfit(G1_Map5temp-G2_Map5temp);
                P_map5 = normcdf(Res1orig.res.ResultMap5-Res2orig.res.ResultMap5,mu,sig);
                DAT = zeros(dims);
                DAT(find(DATMASK)) = P_map5;
                DynamicBC_write_NIFTI(DAT,vmask,Outfilenametemp5)
            else % coef
                for i_ord = 1:RealComp1.RealCompPara.GCAorder
                    
                    Outfilenametemp1 = fullfile(Outputdir,['PermP_x2y_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'.nii']);
                    Outfilenametemp2 = fullfile(Outputdir,['PermP_y2x_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'.nii']);
                  
                    for iperm = 1:Permnum
                        G1_Map1temp(iperm,:) = Res1(iperm).coef.ResultMap1{i_ord};
                        G1_Map2temp(iperm,:) = Res1(iperm).coef.ResultMap2{i_ord};
                        
                        G2_Map1temp(iperm,:) = Res2(iperm).coef.ResultMap1{i_ord};
                        G2_Map2temp(iperm,:) = Res2(iperm).coef.ResultMap2{i_ord};
                    end
                    [mu sig] = normfit(G1_Map1temp-G2_Map1temp);
                    P_map1 = normcdf(Res1orig.coef.ResultMap1{i_ord}-Res2orig.coef.ResultMap1{i_ord},mu,sig);
                    DAT = zeros(dims);
                    DAT(find(DATMASK)) = P_map1;
                    DynamicBC_write_NIFTI(DAT,vmask,Outfilenametemp1)
                    
                    [mu sig] = normfit(G1_Map2temp-G2_Map2temp);
                    P_map2 = normcdf(Res1orig.coef.ResultMap2{i_ord}-Res2orig.coef.ResultMap2{i_ord},mu,sig);
                    DAT = zeros(dims);
                    DAT(find(DATMASK)) = P_map2;
                    DynamicBC_write_NIFTI(DAT,vmask,Outfilenametemp2)
                end
            end
        end
    end
    
elseif strcmp(RealComp1.RealCompPara.mod,'volroi')    
    load(fullfile(Inputdir1,'ROIsignal.mat'));
    ROIsig_g1 = ROIsignals;
    load(fullfile(Inputdir2,'ROIsignal.mat'));
    ROIsig_g2 = ROIsignals;
    
    if Parameter.Permlab==1
        Permnum = Parameter.PermNum;
        G1num = size(ROIsig_g1,1);
        G2num = size(ROIsig_g2,1);
        Totalnum = G1num+G2num;
        totalsig = [ROIsig_g1;ROIsig_g2];
        
        if RealComp1.RealCompPara.COVcond(1)
            COV1 = RealComp1.RealCompPara.COVs;
            COV2 = RealComp2.RealCompPara.COVs;
            covtotal = [COV1;COV2];
        end
        RealComp1.RealCompPara.PermMark = 0;
        RealComp2.RealCompPara.PermMark = 0;
        for i = 1:size(ROIsig_g1,2)
            Res1(i) = BCCT_Perm_CaSCNMap_Cal(RealComp1.RealCompPara,ROIsig_g1,ROIsig_g1(:,i));
            Res2(i) = BCCT_Perm_CaSCNMap_Cal(RealComp2.RealCompPara,ROIsig_g2,ROIsig_g2(:,i));
        end
        for iperm = 1:Permnum
            randord = randperm(Totalnum);
            G1ord = randord(1:G1num);
            G2ord = randord(G1num+1:Totalnum);
            Sig1temp = totalsig(G1ord,:);
            Sig2temp = totalsig(G2ord,:);
            RealComp1temp = RealComp1.RealCompPara;
            RealComp2temp = RealComp2.RealCompPara;
            if RealComp1.RealCompPara.COVcond(1)
                RealComp1temp.COVs = covtotal(G1ord,:);
                RealComp2temp.COVs = covtotal(G2ord,:);
            end
            RealComp1temp.PermMark = 0;
            RealComp2temp.PermMark = 0;
            
            for i = 1:size(Sig1temp,2)
                Res1temp(iperm,i) = BCCT_Perm_CaSCNMap_Cal(RealComp1temp,Sig1temp,Sig1temp(:,i));
                Res2temp(iperm,i) = BCCT_Perm_CaSCNMap_Cal(RealComp2temp,Sig2temp,Sig2temp(:,i));
            end
        end
        if RealComp1.RealCompPara.Calmethod1  % res
            for i = 1:size(Sig1temp,2)
                G1_ResultMap1 = Res1(i).res.ResultMap1;
                G1_ResultMap2 = Res1(i).res.ResultMap2;
                G1_ResultMap3 = Res1(i).res.ResultMap3;
                G1_ResultMap4 = Res1(i).res.ResultMap4;
                G1_ResultMap5 = Res1(i).res.ResultMap5;
                
                G2_ResultMap1 = Res2(i).res.ResultMap1;
                G2_ResultMap2 = Res2(i).res.ResultMap2;
                G2_ResultMap3 = Res2(i).res.ResultMap3;
                G2_ResultMap4 = Res2(i).res.ResultMap4;
                G2_ResultMap5 = Res2(i).res.ResultMap5;
                
                for iperm = 1:Permnum
                    G1_Res1temp1(iperm,:) = Res1temp(iperm,i).res.ResultMap1;
                    G1_Res1temp2(iperm,:) = Res1temp(iperm,i).res.ResultMap2;
                    G1_Res1temp3(iperm,:) = Res1temp(iperm,i).res.ResultMap3;
                    G1_Res1temp4(iperm,:) = Res1temp(iperm,i).res.ResultMap4;
                    G1_Res1temp5(iperm,:) = Res1temp(iperm,i).res.ResultMap5;
                    
                    G2_Res1temp1(iperm,:) = Res2temp(iperm,i).res.ResultMap1;
                    G2_Res1temp2(iperm,:) = Res2temp(iperm,i).res.ResultMap2;
                    G2_Res1temp3(iperm,:) = Res2temp(iperm,i).res.ResultMap3;
                    G2_Res1temp4(iperm,:) = Res2temp(iperm,i).res.ResultMap4;
                    G2_Res1temp5(iperm,:) = Res2temp(iperm,i).res.ResultMap5;
                end
                [mu sig] = normfit(G1_Res1temp1);
                P_map1(i,:) = normcdf(G1_ResultMap1-G2_ResultMap1,mu,sig);
                [mu sig] = normfit(G1_Res1temp2);
                P_map2(i,:) = normcdf(G1_ResultMap2-G2_ResultMap2,mu,sig);
                [mu sig] = normfit(G1_Res1temp3);
                P_map3(i,:) = normcdf(G1_ResultMap3-G2_ResultMap3,mu,sig);
                [mu sig] = normfit(G1_Res1temp4);
                P_map4(i,:) = normcdf(G1_ResultMap4-G2_ResultMap4,mu,sig);
                [mu sig] = normfit(G1_Res1temp5);
                P_map5(i,:) = normcdf(G1_ResultMap5-G2_ResultMap5,mu,sig);
            end
            GroupCompPerm.res.P_mat1 = P_map1;
            GroupCompPerm.res.P_mat2 = P_map2;
            GroupCompPerm.res.P_mat3 = P_map3;
            GroupCompPerm.res.P_mat4 = P_map4;
            GroupCompPerm.res.P_mat5 = P_map5;
            Outnames = fullfile(Outputdir,'PermtestForGCA.mat');
            save(Outnames,'GroupCompPerm');
%             save(Outnames,'P_map1','P_map2','P_map3','P_map4','P_map5');
        else % coef
            for i_ord = 1:RealComp1.RealCompPara.GCAorder
                for i = 1:size(Sig1temp,2)
                    G1_ResultMap1 = Res1(i).coef.ResultMap1{i_ord};
                    G1_ResultMap2 = Res1(i).coef.ResultMap2{i_ord};
                    
                    G2_ResultMap1 = Res2(i).coef.ResultMap1{i_ord};
                    G2_ResultMap2 = Res2(i).coef.ResultMap2{i_ord};
                    
                    for iperm = 1:Permnum
                        G1_Res1temp1(iperm,:) = Res1temp(iperm,i).coef.ResultMap1{i_ord};
                        G1_Res1temp2(iperm,:) = Res1temp(iperm,i).coef.ResultMap2{i_ord};
                        
                        G2_Res1temp1(iperm,:) = Res2temp(iperm,i).coef.ResultMap1{i_ord};
                        G2_Res1temp2(iperm,:) = Res2temp(iperm,i).coef.ResultMap2{i_ord};
                    end
                    [mu sig] = normfit(G1_Res1temp1);
                    P_map1(i,:) = normcdf((G1_ResultMap1-G2_ResultMap1)',mu,sig);
                    [mu sig] = normfit(G1_Res1temp2);
                    P_map2(i,:) = normcdf((G1_ResultMap2-G2_ResultMap2)',mu,sig);
                end
                
                GroupCompPerm.coef.P_mat1 = P_map1;
                GroupCompPerm.coef.P_mat2 = P_map2;
                Outnames = fullfile(Outputdir,['Order_',num2str(i_ord),'_PermtestForGCA.mat']);
                save(Outnames,'GroupCompPerm');
                clear GroupCompPerm
%                 save(Outnames,'P_map1','P_map2');
            end
        end
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
    
    if Parameter.Permlab==1
        Permnum = Parameter.PermNum;
        for i = 1:size(ROIsig_g1,2)    
            Res1orig = BCCT_Perm_CaSCNMap_Cal(RealComp1.RealCompPara,maskedsignals_g1,ROIsig_g1(:,i));
            Res2orig = BCCT_Perm_CaSCNMap_Cal(RealComp2.RealCompPara,maskedsignals_g2,ROIsig_g2(:,i));
            
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
            end
            for iperm = 1:Permnum
                randord = randperm(Totalnum);
                Ord1 = randord(1:G1num);
                Ord2 = randord(G1num+1:Totalnum);
                corrsig1perm = corrsig(Ord1,:);
                corrsig2perm = corrsig(Ord2,:);
                maskedsignal1perm = MASKEDSIGNALS(Ord1,:);
                maskedsignal2perm = MASKEDSIGNALS(Ord2,:);
                RealCompPara1temp = RealComp1.RealCompPara;
                RealCompPara2temp = RealComp2.RealCompPara;
                if RealComp1.RealCompPara.COVcond(1)
                    RealCompPara1temp.COVs = covtotal(Ord1,:);
                    RealCompPara2temp.COVs = covtotal(Ord2,:);
                end
                RealCompPara1temp.PermMark=0;
                RealCompPara2temp.PermMark=0;
                
                Res1(iperm) = BCCT_Perm_CaSCNMap_Cal(RealCompPara1temp,maskedsignal1perm,corrsig1perm(:,i));
                Res2(iperm) = BCCT_Perm_CaSCNMap_Cal(RealCompPara2temp,maskedsignal2perm,corrsig2perm(:,i));

            end
            if RealComp1.RealCompPara.Calmethod1 % res
                
                Outfilenametemp1_L = fullfile(Outputdir,['lh.PermP_x2y_ROI',sprintf('%05d',i),'.mgh']);
                Outfilenametemp2_L = fullfile(Outputdir,['lh.PermP_y2x_ROI',sprintf('%05d',i),'.mgh']);
                Outfilenametemp3_L = fullfile(Outputdir,['lh.PermP_x2ytrans_ROI',sprintf('%05d',i),'.mgh']);
                Outfilenametemp4_L = fullfile(Outputdir,['lh.PermP_y2xtrans_ROI',sprintf('%05d',i),'.mgh']);
                Outfilenametemp5_L = fullfile(Outputdir,['lh.PermP_Netx2y_ROI',sprintf('%05d',i),'.mgh']);
                OutfilenametempMat = fullfile(Outputdir,['PermP_ROI',sprintf('%05d',i),'.mat']);
                
                Outfilenametemp1_R = fullfile(Outputdir,['rh.PermP_x2y_ROI',sprintf('%05d',i),'.mgh']);
                Outfilenametemp2_R = fullfile(Outputdir,['rh.PermP_y2x_ROI',sprintf('%05d',i),'.mgh']);
                Outfilenametemp3_R = fullfile(Outputdir,['rh.PermP_x2ytrans_ROI',sprintf('%05d',i),'.mgh']);
                Outfilenametemp4_R = fullfile(Outputdir,['rh.PermP_y2xtrans_ROI',sprintf('%05d',i),'.mgh']);
                Outfilenametemp5_R = fullfile(Outputdir,['rh.PermP_Netx2y_ROI',sprintf('%05d',i),'.mgh']);
                
                
                for iperm = 1:Permnum
                    G1_Map1temp(iperm,:) = Res1(iperm).res.ResultMap1;
                    G1_Map2temp(iperm,:) = Res1(iperm).res.ResultMap2;
                    G1_Map3temp(iperm,:) = Res1(iperm).res.ResultMap3;
                    G1_Map4temp(iperm,:) = Res1(iperm).res.ResultMap4;
                    G1_Map5temp(iperm,:) = Res1(iperm).res.ResultMap5;
                    G2_Map1temp(iperm,:) = Res2(iperm).res.ResultMap1;
                    G2_Map2temp(iperm,:) = Res2(iperm).res.ResultMap2;
                    G2_Map3temp(iperm,:) = Res2(iperm).res.ResultMap3;
                    G2_Map4temp(iperm,:) = Res2(iperm).res.ResultMap4;
                    G2_Map5temp(iperm,:) = Res2(iperm).res.ResultMap5;
                end
                [mu sig] = normfit(G1_Map1temp-G2_Map1temp);
                P_map1 = normcdf(Res1orig.res.ResultMap1-Res2orig.res.ResultMap1,mu,sig);
                DAT = ones(size(MASK))*0.5;
                DAT(find(MASK)) = P_map1;
                save_mgh(DAT(1:N1size),Outfilenametemp1_L,Ml);
                save_mgh(DAT(1+N1size:N1size+N2size),Outfilenametemp1_R,Mr);
%                 SurfStatWriteData(Outfilenametemp1,DAT)
                
                [mu sig] = normfit(G1_Map2temp-G2_Map2temp);
                P_map2 = normcdf(Res1orig.res.ResultMap2-Res2orig.res.ResultMap2,mu,sig);
                DAT = ones(size(MASK))*0.5;
                DAT(find(MASK)) = P_map2;
                save_mgh(DAT(1:N1size),Outfilenametemp2_L,Ml);
                save_mgh(DAT(1+N1size:N1size+N2size),Outfilenametemp2_R,Mr);
%                 SurfStatWriteData(Outfilenametemp2,DAT)
                
                [mu sig] = normfit(G1_Map3temp-G2_Map3temp);
                P_map3 = normcdf(Res1orig.res.ResultMap3-Res2orig.res.ResultMap3,mu,sig);
                DAT = ones(size(MASK))*0.5;
                DAT(find(MASK)) = P_map3;
                save_mgh(DAT(1:N1size),Outfilenametemp3_L,Ml);
                save_mgh(DAT(1+N1size:N1size+N2size),Outfilenametemp3_R,Mr);
%                 SurfStatWriteData(Outfilenametemp3,DAT)
                
                [mu sig] = normfit(G1_Map4temp-G2_Map4temp);
                P_map4 = normcdf(Res1orig.res.ResultMap4-Res2orig.res.ResultMap4,mu,sig);                
                DAT = ones(size(MASK))*0.5;
                DAT(find(MASK)) = P_map4;
                save_mgh(DAT(1:N1size),Outfilenametemp4_L,Ml);
                save_mgh(DAT(1+N1size:N1size+N2size),Outfilenametemp4_R,Mr);
%                 SurfStatWriteData(Outfilenametemp4,DAT)
                
                [mu sig] = normfit(G1_Map5temp-G2_Map5temp);
                P_map5 = normcdf(Res1orig.res.ResultMap5-Res2orig.res.ResultMap5,mu,sig);              
                DAT = ones(size(MASK))*0.5;
                DAT(find(MASK)) = P_map5;
                save_mgh(DAT(1:N1size),Outfilenametemp5_L,Ml);
                save_mgh(DAT(1+N1size:N1size+N2size),Outfilenametemp5_R,Mr);
%                 SurfStatWriteData(Outfilenametemp5,DAT)
                
                save(OutfilenametempMat,'MASK','P_map1','P_map2','P_map3','P_map4','P_map5')
            else % coef
                for i_ord = 1:RealComp1.RealCompPara.GCAorder
                    
                    Outfilenametemp1_L = fullfile(Outputdir,['lh.PermP_x2y_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'.mgh']);
                    Outfilenametemp2_L = fullfile(Outputdir,['lh.PermP_y2x_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'.mgh']);
                    Outfilenametempmat = fullfile(Outputdir,['PermP_ROI',sprintf('%05d',i),num2str(i_ord),'.mat']);
                    Outfilenametemp1_R = fullfile(Outputdir,['rh.PermP_x2y_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'.mgh']);
                    Outfilenametemp2_R = fullfile(Outputdir,['rh.PermP_y2x_ROI',sprintf('%05d',i),'_Order_',num2str(i_ord),'.mgh']);
                    
                    for iperm = 1:Permnum
                        G1_Map1temp(iperm,:) = Res1(iperm).coef.ResultMap1{i_ord};
                        G1_Map2temp(iperm,:) = Res1(iperm).coef.ResultMap2{i_ord};
                        
                        G2_Map1temp(iperm,:) = Res2(iperm).coef.ResultMap1{i_ord};
                        G2_Map2temp(iperm,:) = Res2(iperm).coef.ResultMap2{i_ord};
                    end
                    [mu sig] = normfit(G1_Map1temp-G2_Map1temp);
                    P_map1 = normcdf(Res1orig.coef.ResultMap1{i_ord}-Res2orig.coef.ResultMap1{i_ord},mu,sig);
                    DAT = ones(size(MASK))*0.5;
                    DAT(find(MASK)) = P_map1;
                    save_mgh(DAT(1:N1size),Outfilenametemp1_L,Ml);
                    save_mgh(DAT(1+N1size:N1size+N2size),Outfilenametemp1_R,Mr);
%                     SurfStatWriteData(Outfilenametemp1,DAT)
                    
                    [mu sig] = normfit(G1_Map2temp-G2_Map2temp);
                    P_map2 = normcdf(Res1orig.coef.ResultMap2{i_ord}-Res2orig.coef.ResultMap2{i_ord},mu,sig);
                    DAT = ones(size(MASK))*0.5;
                    DAT(find(MASK)) = P_map2;
                    save_mgh(DAT(1:N1size),Outfilenametemp2_L,Ml);
                    save_mgh(DAT(1+N1size:N1size+N2size),Outfilenametemp2_R,Mr);
%                     SurfStatWriteData(Outfilenametemp2,DAT)
                    
                    save(Outfilenametempmat,'MASK','P_map1','P_map2')
                end
            end
        end
        
    end
elseif strcmp(RealComp1.RealCompPara.mod,'surfROI')
    disp('comming soon')
end
end