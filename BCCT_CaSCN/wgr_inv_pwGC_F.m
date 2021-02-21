function gc = wgr_inv_pwGC_F(pvalue,nobs,porder);
n = nobs - porder;
th = 1+finv(1-pvalue,porder,n-2*porder-1)/(n-2*porder-1)*porder;
gc = log(th);