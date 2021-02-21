function [x2ycoef,y2xcoef,x2y_T,y2x_T,dfx2y,dfy2x] = gca_coef_sub(AROITimeCourse,ABrain4D,i,theCovariables,Order,nDim4)


xsig = AROITimeCourse;ysig = ABrain4D(:,i);
if std(ysig)~=0
    Y_TimeCourse = ABrain4D(:,i);
    for k = 1:Order,%set order
        AX(:,k)=AROITimeCourse(k:nDim4-Order+k-1);
        BY(:,k)=Y_TimeCourse(k:nDim4-Order+k-1);
    end
    %         Xsig_ux2y = xsig(1:nDim4-Order);
    Xsig_uy2x = xsig(Order+1:nDim4);
    Ysig_ux2y = ysig(Order+1:nDim4);
    %         Ysig_uy2x = ysig(1:nDim4-Order);
    MODAX = term(AX);
    MODBY = term(BY);
    MODX2YC = term([BY,theCovariables]);
    MODY2XC = term([AX,theCovariables]);
    mod1 = 1+MODAX+MODX2YC;
    mod2 = 1+MODBY+MODY2XC;
    mod_x2y = SurfStatLinMod(Ysig_ux2y,mod1);
    mod_y2x = SurfStatLinMod(Xsig_uy2x,mod2);
    if Order>1
        for k = 1:Order
            eval(['mod_x2yt = SurfStatT(mod_x2y,MODAX.AX',num2str(k),');']);
            eval(['mod_y2xt = SurfStatT(mod_y2x,MODBY.BY',num2str(k),');']);
            x2ycoef(1,k) = mod_x2yt.coef(k);
            y2xcoef(1,k) = mod_y2xt.coef(k);
            x2y_T(1,k) = mod_x2yt.t;
            y2x_T(1,k) = mod_y2xt.t;
            x2y_df(1,k) = mod_x2yt.df;
            y2x_df(1,k) = mod_y2xt.df;
        end
    else
        k = 1;
        mod_x2yt = SurfStatT(mod_x2y,MODAX.AX);
        mod_y2xt = SurfStatT(mod_y2x,MODBY.BY);
        x2ycoef(1,k) = mod_x2yt.coef(k);
        y2xcoef(1,k) = mod_y2xt.coef(k);
        x2y_T(1,k) = mod_x2yt.t;
        y2x_T(1,k) = mod_y2xt.t;
        x2y_df(1,k) = mod_x2yt.df;
        y2x_df(1,k) = mod_y2xt.df;
    end
else
    x2ycoef(1,1:Order) = 0;
    y2xcoef(1,1:Order) = 0;
    x2y_T(1,1:Order) = 0;
    y2x_T(1,1:Order) = 0;
    y2x_df = 0;
    x2y_df = 0;
end
df1 = unique(x2y_df);
df2 = unique(y2x_df);
dfx2y = df1(end);
dfy2x = df2(end);
end