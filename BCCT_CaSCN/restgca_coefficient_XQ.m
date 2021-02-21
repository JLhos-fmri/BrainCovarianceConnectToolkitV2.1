function [ResultMap1,ResultMap2,ResultMap3,ResultMap4,res]=restgca_coefficient_XQ(AROITimeCourse, ABrain4D,Order,theCovariables)
         ResultMap1={};ResultMap2={};ResultMap3={};ResultMap4={};
         nDim4 = size(AROITimeCourse,1);
%          [nDim1, nDim2, nDim3, nDim4]=size(ABrain4D); % component by QX
%          AROITimeCourse = reshape(AROITimeCourse, 1, nDim4)';
         AROITimeCourse_now = AROITimeCourse(Order+1:nDim4);
        ndim = size(ABrain4D,2);
%          ABrain4D = reshape(ABrain4D, nDim1*nDim2*nDim3, nDim4)'; % component by QX
         
         theCovariables = [theCovariables(Order+1:end,:),ones(nDim4-Order,1)];

         AX = ones(nDim4-Order,Order);
         BY = ones(nDim4-Order,Order);
   
         AllResultX2Y = ones(Order,ndim);AllResultX2Y_AR = ones(Order,ndim);
         AllResultY2X = ones(Order,ndim);AllResultY2X_AR = ones(Order,ndim);
         
         for i = 1:ndim,
             Y_TimeCourse = ABrain4D(:,i);
             Y_TimeCourse_now = Y_TimeCourse(Order+1:nDim4);
             theJudgment = var(Y_TimeCourse_now);
             
             for k = 1:Order,%set order       
                 AX(:,k)=AROITimeCourse(k:nDim4-Order+k-1);
                 BY(:,k)=Y_TimeCourse(k:nDim4-Order+k-1);
             end
             meanAX = sqrt(sum(AX.^2));
             stdAX = std(AX);
             meanBY = sqrt(sum(BY.^2));
             stdBY = std(BY);
             Regressors = [AX,BY,theCovariables];
             ResultX2Y = rest_regress(Y_TimeCourse_now,Regressors);
             ResultY2X = rest_regress(AROITimeCourse_now,Regressors);
             AllResultX2Y(:,i) = ResultX2Y(1:Order)';
             AllResultY2X(:,i) = ResultY2X(Order+1:Order*2)';
             AllResultX2Y_AR(:,i) = ResultX2Y(Order+1:Order*2)';
             
             for iord = 1:Order
                 T_X2Y(iord) = ResultX2Y(iord)*meanAX(iord)/stdAX(iord);
                 T_Y2X(iord) = ResultY2X(iord)*meanBY(iord)/stdBY(iord);
             end
             AllResultX2YT(:,i) = T_X2Y;
             AllResultY2XT(:,i) = T_Y2X;
%              size(Regressors,1)
             [Z1(:,i), P1(:,i)] = AS_TFRtoZ(T_X2Y,'T',size(Regressors,1)-2,[]); 
             [Z2(:,i), P2(:,i)] = AS_TFRtoZ(T_Y2X,'T',size(Regressors,1)-2,[]); 
             
             if theJudgment == 0,
                AllResultY2X_AR(:,i) = zeros(Order,1);
             else
                AllResultY2X_AR(:,i) = ResultY2X(1:Order)';
             end
             
         end
         
         for j = 1:Order
             res.T_T1{j} = AllResultX2YT(j,:);
             res.T_T2{j} = AllResultY2XT(j,:);
             res.Z_T1{j} = Z1(j,:);
             res.Z_T2{j} = Z2(j,:);
             res.P_T1{j} = P1(j,:);
             res.P_T2{j} = P2(j,:);
         end
         
         for j = 1:Order,
             ResultMap1 = [ResultMap1,{AllResultX2Y(j,:)}];
             ResultMap2 = [ResultMap2,{AllResultY2X(j,:)}];
             ResultMap3 = [ResultMap3,{AllResultX2Y_AR(j,:)}];
             ResultMap4 = [ResultMap4,{AllResultY2X_AR(j,:)}];       
         end
end

function beta = rest_regress(y,X)
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
end