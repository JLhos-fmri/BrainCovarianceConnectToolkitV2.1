function BCCT_MatrixGCshow_coef_GUI(RESshow,outdir,enhanind,Coltype,PVALT,Pvalused,orders)

% load(indir);
[mpat,mnam,mext] = fileparts(which('BCCT_MatrixGCshow_coef_GUI.m'));
load(fullfile(mpat,'gray.mat'));
load(fullfile(mpat,'MatColormap','mycolormap.mat'));

GCAorder = size(RESshow,2);
for iorder = 1:GCAorder
    X2Y = RESshow(iorder).GCA_coef_x2y;
    Y2X = RESshow(iorder).GCA_coef_y2x;
    X2Y_T = RESshow(iorder).GCA_coef_x2yT;
    Y2X_T = RESshow(iorder).GCA_coef_y2xT;
    X2Y_Z = RESshow(iorder).GCA_coef_x2yZ;
    Y2X_Z = RESshow(iorder).GCA_coef_y2xZ;
    P_X2Y = RESshow(iorder).GCA_coef_x2y_pval;
    P_Y2X = RESshow(iorder).GCA_coef_y2x_pval;
    PermLab = 0;
    if isfield(RESshow(1),'GCA_coef_x2y_permpval')
        PermLab = 1;
        Perm_X2Y = RESshow.GCA_coef_x2y_permpval;
        Perm_Y2X = RESshow.GCA_coef_y2x_permpval;
    end
    %
    X2Y(1:size(X2Y,1)+1:end) = 0;
    X2Y = X2Y(orders,orders);
    Y2X(1:size(Y2X,1)+1:end) = 0;
    Y2X = Y2X(orders,orders);
    
    X2Y_T(1:size(X2Y_T,1)+1:end) = 0;
    X2Y_T = X2Y_T(orders,orders);
    Y2X_T(1:size(Y2X_T,1)+1:end) = 0;
    Y2X_T = Y2X_T(orders,orders);
    
    X2Y_Z(1:size(X2Y_Z,1)+1:end) = 0;
    X2Y_Z = X2Y_Z(orders,orders);
    Y2X_Z(1:size(Y2X_Z,1)+1:end) = 0;
    Y2X_Z = Y2X_Z(orders,orders);
    
    P_X2Y(1:size(P_X2Y,1)+1:end) = 0;
    P_X2Y = P_X2Y(orders,orders);
    P_Y2X(1:size(P_Y2X,1)+1:end) = 0;
    P_Y2X = P_Y2X(orders,orders);
  
    if PermLab
        Perm_X2Y(1:size(Perm_X2Y,1)+1:end) = 0;
        Perm_X2Y = Perm_X2Y(orders,orders);
        Perm_Y2X(1:size(Perm_Y2X,1)+1:end) = 0;
        Perm_Y2X = Perm_Y2X(orders,orders);
    end
    
    %
    NROI = size(X2Y,1);
    if NROI>50
        width1 = 0.01;
        width2 = 1;
    else
        width1 = 0.5;
        width2 = 3;
    end
    uesdmat = ones(NROI);
    uesdmat(1:NROI+1:end) = 0;
    %% set Pvalues
    if PVALT==1
        PVALSHOW = Pvalused;
    elseif PVALT==2
        PVALSHOW = 1/(NROI*(NROI-1)/2);
    elseif PVALT==3
        Pusedinfo = P_X2Y(find(uesdmat));
        [h pi] = fdr(Pusedinfo, Pvalused);
        hind = find(h);
        if ~isempty(hind)
            PVALSHOW = max(pi(hind));
        else
            PVALSHOW = 1e-30;
        end
        if PermLab
            Pusedinfo = 1-Perm_X2Y(find(uesdmat));
            [h pi] = fdr(Pusedinfo, Pvalused);
            hind = find(h);
            if ~isempty(hind)
                PVALSHOW_PermX2Y = max(pi(hind));
            else
                PVALSHOW_PermX2Y = 1e-30;
            end
            
            Pusedinfo = 1-Perm_Y2X(find(uesdmat));
            [h pi] = fdr(Pusedinfo, Pvalused);
            hind = find(h);
            if ~isempty(hind)
                PVALSHOW_PermY2X = max(pi(hind));
            else
                PVALSHOW_PermY2X = 1e-30;
            end
        end
    else
        PVALSHOW = Pvalused/(NROI*(NROI-1)/2);
    end
    if PVALSHOW<1e-6
        PVALSHOW1 = 1e-6;
    else
        PVALSHOW1 = PVALSHOW;
    end
    PSHOW = find((P_X2Y)<=PVALSHOW*2);
    P3 = zeros(size(P_X2Y));
    P3(PSHOW) = 1;
    
    colormapshow0 = AFNICOLORMAP(128);
    if Coltype==1
        colormapshow(1:64,:) = colormapshow0(1:64,:);
        colormapshow(66:129,:) = colormapshow0(65:end,:);
        colormapshow(65,:) = [1 1 1];
    elseif Coltype==2
        colormapshow = RdBu_map(end:-1:1,:);
        colormapshow(128,:) = [1 1 1];
        colormapshow(129,:) = [1 1 1];
    else
        colormapshow = RdYlBu_map(end:-1:1,:);
        colormapshow(128,:) = [1 1 1];
        colormapshow(129,:) = [1 1 1];
    end
    enhanind = unique([0,NROI,enhanind]);
    
    %% Figure Show
    % save test
    Show_GC_coef_subfun(X2Y,P3,NROI,colormapshow,enhanind,outdir,['Order-',num2str(iorder),'-CaSCNcoefX2Y'],width1,width2);
    Show_GC_coef_subfun(Y2X,P3',NROI,colormapshow,enhanind,outdir,['Order-',num2str(iorder),'-CaSCNcoefY2X'],width1,width2);
    
    Show_GC_coef_subfun(X2Y_T,P3,NROI,colormapshow,enhanind,outdir,['Order-',num2str(iorder),'-CaSCNcoefX2Y-T'],width1,width2);
    Show_GC_coef_subfun(Y2X_T,P3',NROI,colormapshow,enhanind,outdir,['Order-',num2str(iorder),'-CaSCNcoefY2X-T'],width1,width2);
    
    Show_GC_coef_subfun(X2Y_Z,P3,NROI,colormapshow,enhanind,outdir,['Order-',num2str(iorder),'-CaSCNcoefX2Y-Z'],width1,width2);
    Show_GC_coef_subfun(Y2X_Z,P3',NROI,colormapshow,enhanind,outdir,['Order-',num2str(iorder),'-CaSCNcoefY2X-Z'],width1,width2);
    
    if PermLab
        if PVALT==3
            PSHOW = find((1-Perm_X2Y)<=PVALSHOW_PermX2Y|(Perm_X2Y<=PVALSHOW_PermX2Y&Perm_X2Y>0));
        else
            PSHOW = find((1-Perm_X2Y)<=PVALSHOW|(Perm_X2Y<=PVALSHOW&Perm_X2Y>0));
        end
        P3 = zeros(size(Perm_X2Y));
        P3(PSHOW) = 1;
        Show_GC_coef_subfun(X2Y,P3,NROI,colormapshow,enhanind,outdir,['Order-',num2str(iorder),'-CaSCNcoefX2Y_Perm'],width1,width2);
        
        if PVALT==3
            PSHOW = find((1-Perm_Y2X)<=PVALSHOW_PermY2X|(Perm_Y2X<=PVALSHOW_PermY2X&Perm_Y2X>0));
        else
            PSHOW = find((1-Perm_Y2X)<=PVALSHOW|(Perm_Y2X<=PVALSHOW&Perm_Y2X>0));
        end
        P3 = zeros(size(Perm_Y2X));
        P3(PSHOW) = 1;
        Show_GC_coef_subfun(Y2X,P3,NROI,colormapshow,enhanind,outdir,['Order-',num2str(iorder),'-CaSCNcoefY2X_Perm'],width1,width2);
        
    end
end
end