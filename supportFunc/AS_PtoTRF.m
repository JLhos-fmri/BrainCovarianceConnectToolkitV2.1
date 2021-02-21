function [Z,out] = AS_PtoTRF(pval,Flag,Df1,Df2)
% This is from software rest.
% The orginal version is [Z P Header] = rest_TFRtoZ(ImgFile,OutputName,Flag,Df1,Df2,Header)

% The bellow is the orginal annotation.
% function [Z P Header] = rest_TFRtoZ(ImgFile,OutputName,Flag,Df1,Df2,Header)
% FORMAT [Z P Header] = y_TFRtoZ(ImgFile,OutputName,Flag,Df1,Df2,Header)
%   Input:
%     ImgFile    - T, F or R statistical image which wanted to be converted to Z statistical value
%     OutputName - The output name
%     Flag       - 'T', 'F' or 'R'. Indicate the type of the input statsical image
%                - If not defined or defined as empty, then will read the statistical type and degree of freedom information from the image (if the statistical analysis was performed with REST or SPM).
%     Df1        - the degree of freedom of the statistical image. For F statistical image, there is also Df2
%     Df2        - the second degree of freedom of F statistical image
%     Header     - If ImgFile is a MATRIX, the Nifti Header should be specified 
%   Output:
%     Z          - Z statistical image. Also output as .img/.hdr.
%     P          - The corresponding P value
%     Header     - the output Nifti Header
%___________________________________________________________________________
% Written by YAN Chao-Gan 100424.
% State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University, China, 100875
% ycg.yan@gmail.com
% Revised by YAN Chao-Gan 100814. Changed to call spm_t2z for T and R images to use approximation in case of big T values.
% Last revised by YAN Chao-Gan 121223. Convert the big F values to Z values in an approximation like spm_t2z to treat big T values. Detect the Flag and DF from the data if Flag is not defined.

% For me, I use the idea and rewrite some of it.

%Added by YAN Chao-Gan 121222. Detect the Flag and DF from the data if Flag is not defined.


if strcmpi(Flag,'F')
    fprintf('Convert F to Z...\n');
    Data = 0.01:0.00001:100;
    P = 1-fcdf(Data,Df1,Df2);
    ind = find(P<pval);
    out = Data(ind(1)-1);
    Z = norminv(fcdf(out,Df1,Df2)); %YAN Chao-Gan 100814. Use one-tail because F value is positive and one-tail.
%     Z(Data==0) = 0;
    

    %YAN Chao-Gan, 121223. Convert the big F values to Z values in an approximation like spm_t2z to treat big T values.
    %Referenced from spm_t2z.m
    
    %Tol = 1E-16; %minimum tail probability for direct computation
    Tol = 1E-10; %minimum tail probability for direct computation. This is the tolorance value used in spm_t2z.
    
    F1    = finv(1 - Tol,Df1,Df2);
    %mQb   = Data > F1;
    mQb   = isinf(Z); %Only deal with those with Inf values. YAN Chao-Gan, 121223.
    if any(mQb(:))
        z1          = -norminv(Tol);
        F2          = F1-[1:5]/10;
        z2          = norminv(fcdf(F2,Df1,Df2));
        %-least squares line through ([f1,t2],[z1,z2]) : z = m*f + c
        mc          = [[F1,F2]',ones(length([F1,F2]),1)] \ [z1,z2]';
        
        %-------------------------------------------------------------------
        %-Logarithmic extrapolation
        %-------------------------------------------------------------------
        l0=1/mc(1);
        %-Perform logarithmic extrapolation, negate z for positive t-values
        Z(mQb) = ( log( Data(mQb) -F1 + l0 ) + (z1-log(l0)) );
        %-------------------------------------------------------------------
    end
    
    
else % T image or R image: YAN Chao-Gan 100814. Changed to call spm_t2z to use approximation in case of big T values.
    
    if strcmpi(Flag,'R')
        Data0 = 0.01:0.00001:1;
%         fprintf('Convert R to T...\n');
        Data = Data0 .* sqrt(Df1 ./ (1 - Data0.*Data0));
    else
        Data = 0.01:0.00001:100;
    end
    
%     fprintf('Convert T to Z...\n');
    
    P = 2*(1-tcdf(abs(Data),Df1)); %Two-tailed P value
    
    ind = find(P<pval);
    if strcmpi(Flag,'R')
        out = Data0(ind(1)-1);
    else
        out = Data(ind(1)-1);
    end
    Tol = 1E-16; %minimum tail probability for direct computation
    Z = spm_t2z(Data(ind(1)-1),Df1,Tol);
end


Z(isnan(Z))=0;
P(isnan(P))=1;

% fprintf('\n\tT/F/R to Z Calculation finished.\n');