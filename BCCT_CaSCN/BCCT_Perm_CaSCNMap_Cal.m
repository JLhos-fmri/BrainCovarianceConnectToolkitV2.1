function Res = BCCT_Perm_CaSCNMap_Cal(RealCompPara,maskedSignal,ROIsignals1)
% ROIsignals1: 1*N
Calmethod1 = RealCompPara.Calmethod1;
COVcond = RealCompPara.COVcond;
if COVcond(1)
    RealCompPara.COVs = zscore(RealCompPara.COVs);
% else
%     RealCompPara.COVs = ones(size(ROIsignals1));
end
ROIsignals1 = zscore(ROIsignals1);
maskedSignal = zscore(maskedSignal);
if COVcond(1)==1
    theCovariables = [ones(size(ROIsignals1,1),1),RealCompPara.COVs];
else
    theCovariables = ones(size(ROIsignals1,1),1);
end
[ix,IOrds] = sort(RealCompPara.IMGORD);
Order = RealCompPara.GCAorder;
if Calmethod1 % res
    [ResultMap1,ResultMap2,ResultMap3,ResultMap4,ResultMap5] = restgca_residual(ROIsignals1(IOrds,:), maskedSignal(IOrds,:),Order,theCovariables(IOrds,:));
    ResultMap_p1 = wgr_pwGC_F(ResultMap1,size(ROIsignals1,1),Order);
    ResultMap_p2 = wgr_pwGC_F(ResultMap2,size(ROIsignals1,1),Order);
    Res.res.ResultMap1 = ResultMap1;
    Res.res.ResultMap2 = ResultMap2;
    Res.res.ResultMap3 = ResultMap3;
    Res.res.ResultMap4 = ResultMap4;
    Res.res.ResultMap5 = ResultMap5;
    Res.res.ResultMap_p1 = ResultMap_p1;
    Res.res.ResultMap_p2 = ResultMap_p2;
else % coef
%     [ResultMap1,ResultMap2,res]=gca_coef_new(AROITimeCourse, ABrain4D,Order,theCovariables)
    [ResultMap1,ResultMap2,res]=gca_coef_new(ROIsignals1(IOrds,:), maskedSignal(IOrds,:),Order,theCovariables(IOrds,:));
%     [ResultMap1,ResultMap2,ResultMap3,ResultMap4,res]=restgca_coefficient_XQ(ROIsignals1(IOrds,:), maskedSignal(IOrds,:),Order,theCovariables(IOrds,:));
    Res.coef.ResultMap1 = ResultMap1;
    Res.coef.ResultMap2 = ResultMap2;
%     Res.coef.ResultMap3 = ResultMap3;
%     Res.coef.ResultMap4 = ResultMap4;
    Res.coef.res = res;
end

if RealCompPara.PermMark
    PermNum = RealCompPara.PermNum;
    
    totallen = size(maskedSignal,2);
    piece = 1000;
    Nlen = ceil(totallen/piece);
%     DataOut1 = zeros(size(DATMASK));
%     DataOut2 = zeros(size(DATMASK));
%     DataOut3 = zeros(size(DATMASK));
%     DataOut4 = zeros(size(DATMASK));
%     DataOut5 = zeros(size(DATMASK));
    for ilen = 1:Nlen
        if ilen~=Nlen
            PIECEORD = 1+piece*(ilen-1):piece*ilen;
        else
            PIECEORD = 1+piece*(ilen-1):totallen;
        end
        clear ResultMap1 ResultMap2 ResultMap3 ResultMap4 ResultMap5 res
        
        for iperm = 1:PermNum
            IOrds = randperm(size(ROIsignals1,1));
            if Calmethod1 % res
                [ResultMap1(iperm,:),ResultMap2(iperm,:),ResultMap3(iperm,:),ResultMap4(iperm,:),ResultMap5(iperm,:)] = restgca_residual(ROIsignals1(IOrds,:), maskedSignal(IOrds,PIECEORD),Order,theCovariables(IOrds,:));
            else
                [ResultMap1(iperm,:),ResultMap2(iperm,:),rest]=gca_coef_new(ROIsignals1(IOrds,:), maskedSignal(IOrds,:),Order,theCovariables(IOrds,:));
%                 [ResultMap1(iperm,:),ResultMap2(iperm,:),ResultMap3(iperm,:),ResultMap4(iperm,:),~]=restgca_coefficient_XQ(ROIsignals1(IOrds,:), maskedSignal(IOrds,PIECEORD),Order,theCovariables(IOrds,:));
            end
        end
        if Calmethod1
            [mu,sig,~,sigci] = normfit(ResultMap1);
            P_map1(PIECEORD) = normcdf(Res.res.ResultMap1(PIECEORD),mu,sig);
            
            [mu,sig,~,sigci] = normfit(ResultMap2);
            P_map2(PIECEORD) = normcdf(Res.res.ResultMap2(PIECEORD),mu,sig);
            
            [mu,sig,~,sigci] = normfit(ResultMap3);
            P_map3(PIECEORD) = normcdf(Res.res.ResultMap3(PIECEORD),mu,sig);
            
            [mu,sig,~,sigci] = normfit(ResultMap4);
            P_map4(PIECEORD) = normcdf(Res.res.ResultMap4(PIECEORD),mu,sig);
            
            [mu,sig,~,sigci] = normfit(ResultMap5);
            P_map5(PIECEORD) = normcdf(Res.res.ResultMap5(PIECEORD),mu,sig);
        else
            for iord = 1:Order
                for iperm = 1:PermNum
                    ResultMap1t(iperm,:) = ResultMap1{iperm,iord};
                    ResultMap2t(iperm,:) = ResultMap2{iperm,iord};
                end
                
                [mu,sig,~,sigci] = normfit(ResultMap1t);
                P_map1(iord,PIECEORD) = normcdf(Res.coef.ResultMap1{iord}(PIECEORD)',mu,sig);
                [mu,sig,~,sigci] = normfit(ResultMap2t);
                P_map2(iord,PIECEORD) = normcdf(Res.coef.ResultMap2{iord}(PIECEORD)',mu,sig);
            end
            clear ResultMap1t ResultMap2t
        end
    end
    if Calmethod1
        Res.res.P_map1 = P_map1;
        Res.res.P_map2 = P_map2;
        Res.res.P_map3 = P_map3;
        Res.res.P_map4 = P_map4;
        Res.res.P_map5 = P_map5;
    else
        Res.coef.P_map1 = P_map1;
        Res.coef.P_map2 = P_map2;
    end
end
end
