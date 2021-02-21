%%%% Compute the GCA on residual %%%%%%%
function [ResultMap1,ResultMap2,ResultMap3,ResultMap4,ResultMap5] = restgca_residual(AROITimeCourse, ABrain4D,Order,theCovariables)
%          [nDim1, nDim2, nDim3, nDim4]=size(ABrain4D); % component by QX
        nDim4 = size(AROITimeCourse,1);
%          AROITimeCourse = reshape(AROITimeCourse, 1, nDim4)';
         AROITimeCourse_now = AROITimeCourse(Order+1:nDim4);
   
%          ABrain4D = reshape(ABrain4D, nDim1*nDim2*nDim3, nDim4)'; % component by QX
         ndim = size(ABrain4D,2);
         theCovariables = [theCovariables(Order+1:end,:),ones(nDim4-Order,1)];
   
         AX = ones(nDim4-Order,Order);
         BY = ones(nDim4-Order,Order);

         AllResult1 = ones(1,ndim);
         AllResult2 = ones(1,ndim);
   
         for i = 1:ndim,
             Y_TimeCourse = ABrain4D(:,i);
             Y_TimeCourse_now = Y_TimeCourse(Order+1:nDim4);
             Judgment = var(Y_TimeCourse_now);
        
             for k = 1:Order,%set order       
                 AX(:,k) = AROITimeCourse(k:nDim4-Order+k-1);
                 BY(:,k) = Y_TimeCourse(k:nDim4-Order+k-1);
             end
        
             Regressors1 = [BY,theCovariables];
             Regressors2 = [AX,theCovariables];
             Regressors3 = [AX,BY,theCovariables];
        
             Residual_Y = rest_regress(Y_TimeCourse_now,Regressors1);
             Residual_X = rest_regress(AROITimeCourse_now,Regressors2);
             Residual_X2Y = rest_regress(Y_TimeCourse_now,Regressors3);
             Residual_Y2X = rest_regress(AROITimeCourse_now,Regressors3);
        
             if Judgment == 0,
                F_X2Y = 0;
                F_Y2X = 0;
             else
                F_X2Y = log(Residual_Y/Residual_X2Y);
                F_Y2X = log(Residual_X/Residual_Y2X);
             end
        
             AllResult1(:,i) = F_X2Y;
             AllResult2(:,i) = F_Y2X;       
         end
    
         ResultMap1 = AllResult1;
         ResultMap2 = AllResult2;
         ResultMap3 = abs((ResultMap1.*(nDim4-Order)-(Order-1)/3)).^0.5;
         ResultMap4 = abs((ResultMap2.*(nDim4-Order)-(Order-1)/3)).^0.5;
         ResultMap5 = ResultMap1-ResultMap2;
end
function Residual_var = rest_regress(y,X)
         [n,ncolX] = size(X);
         [Q,R,perm] = qr(X,0);
         p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));       
         if p < ncolX,
            R = R(1:p,1:p);
            Q = Q(:,1:p);
            perm = perm(1:p);
         end
         beta = zeros(ncolX,1);
         beta(perm) = R \ (Q'*y);
% compute the var of residual           
         yhat = X*beta;                     
         residual = y-yhat;
         Residual_var = sum(residual.^2)/(length(residual)-1);
end