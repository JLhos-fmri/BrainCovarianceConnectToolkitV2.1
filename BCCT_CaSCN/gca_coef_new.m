function [ResultMap1,ResultMap2,res]=gca_coef_new(AROITimeCourse, ABrain4D,Order,theCovariables)
% [ResultMap1,ResultMap2,res]=gca_coef_new(ROIsignals1(IOrds,:), maskedSignal(IOrds,:),Order,theCovariables(IOrds,:))
ResultMap1={};ResultMap2={};
nDim4 = size(AROITimeCourse,1);
AROITimeCourse_now = AROITimeCourse(Order+1:nDim4);
ndim = size(ABrain4D,2);

theCovariables = [theCovariables(Order+1:end,:),ones(nDim4-Order,1)];
AX = ones(nDim4-Order,Order);
BY = ones(nDim4-Order,Order);
% MODCOV = term(theCovariables);
parfor i = 1:ndim
    [x2ycoeft,y2xcoeft,x2y_Tt,y2x_Tt,dfx2yt,dfy2xt] = gca_coef_sub(AROITimeCourse,ABrain4D,i,theCovariables,Order,nDim4);
    x2ycoef(i,:) = x2ycoeft;
    y2xcoef(i,:) = y2xcoeft;
    x2y_T(i,:) = x2y_Tt;
    y2x_T(i,:) = y2x_Tt;
    dfx2y(i,:) = dfx2yt;
    dfy2x(i,:) = dfy2xt;
end
x2y_df = unique(dfx2y);
y2x_df = unique(dfy2x);
for j = 1:length(x2y_df)
    nums(j) = nnz(dfx2y==x2y_df(j));    
end
[maxn, maxind] = max(nums);
[Z1,P1] = AS_TFRtoZ(x2y_T,'T',x2y_df(maxind),[]);

for j = 1:length(y2x_df)
    nums(j) = nnz(dfy2x==y2x_df(j));    
end
[maxn, maxind] = max(nums);
[Z2,P2] = AS_TFRtoZ(y2x_T,'T',y2x_df(maxind),[]);
% for i = 1:ndim,
%     xsig = AROITimeCourse;ysig = ABrain4D(:,i);
%     if std(ysig)~=0
%         Y_TimeCourse = ABrain4D(:,i);
%         for k = 1:Order,%set order
%             AX(:,k)=AROITimeCourse(k:nDim4-Order+k-1);
%             BY(:,k)=Y_TimeCourse(k:nDim4-Order+k-1);
%         end
% %         Xsig_ux2y = xsig(1:nDim4-Order);
%         Xsig_uy2x = xsig(Order+1:nDim4);
%         Ysig_ux2y = ysig(Order+1:nDim4);
% %         Ysig_uy2x = ysig(1:nDim4-Order);
%         MODAX = term(AX);
%         MODBY = term(BY);
%         MODX2YC = term([BY,theCovariables]);
%         MODY2XC = term([AX,theCovariables]);
%         mod1 = 1+MODAX+MODX2YC;
%         mod2 = 1+MODBY+MODY2XC;
%         mod_x2y = SurfStatLinMod(Ysig_ux2y,mod1);
%         mod_y2x = SurfStatLinMod(Xsig_uy2x,mod2);
%         
%         for k = 1:Order
%             eval(['mod_x2yt = SurfStatT(mod_x2y,MODAX.AX',num2str(k),');']);
%             eval(['mod_y2xt = SurfStatT(mod_y2x,MODBY.BY',num2str(k),');']);
%             x2ycoef(i,k) = mod_x2yt.coef(k);
%             y2xcoef(i,k) = mod_y2xt.coef(k);
%             x2y_T(i,k) = mod_x2yt.t;
%             y2x_T(i,k) = mod_y2xt.t;
%             x2y_df(i,k) = mod_x2yt.df;
%             y2x_df(i,k) = mod_y2xt.df;
%         end
%     else
%         x2ycoef(i,1:Order) = 0;
%         y2xcoef(i,1:Order) = 0;
%         x2y_T(i,1:Order) = 0;
%         y2x_T(i,1:Order) = 0;
%     end
% end

for j = 1:Order
    res.T_T1{j} = x2y_T(:,j);
    res.T_T2{j} = y2x_T(:,j);
    res.Z_T1{j} = Z1(:,j);
    res.Z_T2{j} = Z2(:,j);
    res.P_T1{j} = P1(:,j);
    res.P_T2{j} = P2(:,j);
end

for j = 1:Order,
    ResultMap1 = [ResultMap1,{x2ycoef(:,j)}];
    ResultMap2 = [ResultMap2,{y2xcoef(:,j)}];
end
end