function pvalue = wgr_pwGC_F(gc,nobs,porder)
n = nobs - porder;
TH = exp(gc);
THnew = TH-1;
THnew2 = THnew/porder*(n-2*porder-1);
pvalue = fcdf(THnew2,porder,n-2*porder-1);
end
% th = 1+finv(1-pvalue,porder,n-2*porder-1)/(n-2*porder-1)*porder;
% gc = log(th);