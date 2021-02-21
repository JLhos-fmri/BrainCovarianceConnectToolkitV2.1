function BCCT_MatrixGCshow_res_GUI(RESshow,outdir,enhanind,Coltype,PVALT,Pvalused,orders)

% load(indir);
[mpat,mnam,mext] = fileparts(which('BCCT_MatrixGCshow_res_GUI.m'));
load(fullfile(mpat,'gray.mat'));
load(fullfile(mpat,'MatColormap','mycolormap.mat'));

X2Y = RESshow.GCA_res_x2y;
Y2X = RESshow.GCA_res_y2x;
X2Y_trans = RESshow.GCA_res_x2y_trans;
Y2X_trans = RESshow.GCA_res_y2x_trans;
P_X2Y = RESshow.GCA_res_x2y_pval;
P_Y2X = RESshow.GCA_res_y2x_pval;
FX2Y = RESshow.GCA_res_fx2y;
PermLab = 0;
if isfield(RESshow,'GCA_res_x2y_permpval')
    PermLab = 1;
    Perm_X2Y = RESshow.GCA_res_x2y_permpval;
    Perm_Y2X = RESshow.GCA_res_y2x_permpval;
    Perm_X2Y_trans = RESshow.GCA_res_x2y_trans_permpval;
    Perm_Y2X_trans = RESshow.GCA_res_y2x_trans_permpval;
    Perm_FX2Y = RESshow.GCA_res_fx2y_permpval;
end
%
X2Y(1:size(X2Y,1)+1:end) = 0;
X2Y = X2Y(orders,orders);
Y2X(1:size(Y2X,1)+1:end) = 0;
Y2X = Y2X(orders,orders);
X2Y_trans(1:size(X2Y_trans,1)+1:end) = 0;
X2Y_trans = X2Y_trans(orders,orders);
Y2X_trans(1:size(Y2X_trans,1)+1:end) = 0;
Y2X_trans = Y2X_trans(orders,orders);
P_X2Y(1:size(P_X2Y,1)+1:end) = 0;
P_X2Y = P_X2Y(orders,orders);
P_Y2X(1:size(P_Y2X,1)+1:end) = 0;
P_Y2X = P_Y2X(orders,orders);
FX2Y(1:size(FX2Y,1)+1:end) = 0;
FX2Y = FX2Y(orders,orders);
if PermLab    
    Perm_X2Y(1:size(Perm_X2Y,1)+1:end) = 0;
    Perm_X2Y = Perm_X2Y(orders,orders);
    Perm_Y2X(1:size(Perm_Y2X,1)+1:end) = 0;
    Perm_Y2X = Perm_Y2X(orders,orders);
    Perm_X2Y_trans(1:size(Perm_X2Y_trans,1)+1:end) = 0;
    Perm_X2Y_trans = Perm_X2Y_trans(orders,orders);
    Perm_Y2X_trans(1:size(Perm_Y2X_trans,1)+1:end) = 0;
    Perm_Y2X_trans = Perm_Y2X_trans(orders,orders);
    Perm_FX2Y(1:size(Perm_FX2Y,1)+1:end) = 0;
    Perm_FX2Y = Perm_FX2Y(orders,orders);
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
    Pusedinfo = 1-P_X2Y(find(uesdmat));
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
        
        Pusedinfo = 1-Perm_X2Y_trans(find(uesdmat));
        [h pi] = fdr(Pusedinfo, Pvalused);
        hind = find(h);
        if ~isempty(hind)
            PVALSHOW_PermX2Ytrans = max(pi(hind));
        else
            PVALSHOW_PermX2Ytrans = 1e-30;
        end
        
        Pusedinfo = 1-Perm_Y2X_trans(find(uesdmat));
        [h pi] = fdr(Pusedinfo, Pvalused);
        hind = find(h);
        if ~isempty(hind)
            PVALSHOW_PermY2Xtrans = max(pi(hind));
        else
            PVALSHOW_PermY2Xtrans = 1e-30;
        end
        
        Pusedinfo = 1-Perm_FX2Y(find(uesdmat));
        [h pi] = fdr(Pusedinfo, Pvalused);
        hind = find(h);
        if ~isempty(hind)
            PVALSHOW_PermFX2Y = max(pi(hind));
        else
            PVALSHOW_PermFX2Y = 1e-30;
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
PSHOW = find((1-P_X2Y)<=PVALSHOW);
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
Show_GC_res_subfun(X2Y,P3,NROI,colormapshow,enhanind,outdir,'CaSCNresX2Y',width1,width2);
Show_GC_res_subfun(Y2X,P3',NROI,colormapshow,enhanind,outdir,'CaSCNresY2X',width1,width2);

Show_GC_res_subfun(X2Y_trans,P3,NROI,colormapshow,enhanind,outdir,'CaSCNresX2Y_trans',width1,width2);
Show_GC_res_subfun(Y2X_trans,P3',NROI,colormapshow,enhanind,outdir,'CaSCNresY2X_trans',width1,width2);

if PermLab
    if PVALT==3
        PSHOW = find((1-Perm_X2Y)<=PVALSHOW_PermX2Y);
    else
        PSHOW = find((1-Perm_X2Y)<=PVALSHOW);
    end
    P3 = zeros(size(Perm_X2Y));
    P3(PSHOW) = 1;
    Show_GC_res_subfun(X2Y,P3,NROI,colormapshow,enhanind,outdir,'CaSCNresX2Y_Perm',width1,width2);
    
    if PVALT==3
        PSHOW = find((1-Perm_Y2X)<=PVALSHOW_PermY2X);
    else
        PSHOW = find((1-Perm_Y2X)<=PVALSHOW);
    end    
    P3 = zeros(size(Perm_Y2X));
    P3(PSHOW) = 1;
    Show_GC_res_subfun(Y2X,P3,NROI,colormapshow,enhanind,outdir,'CaSCNresY2X_Perm',width1,width2);
    if PVALT==3
        PSHOW = find((1-Perm_X2Y_trans)<=PVALSHOW_PermX2Ytrans);
    else
        PSHOW = find((1-Perm_X2Y_trans)<=PVALSHOW);
    end
    P3 = zeros(size(Perm_X2Y_trans));
    P3(PSHOW) = 1;
    Show_GC_res_subfun(X2Y_trans,P3,NROI,colormapshow,enhanind,outdir,'CaSCNresX2Y_trans_perm',width1,width2);
    
    if PVALT==3
        PSHOW = find((1-Perm_Y2X_trans)<=PVALSHOW_PermY2Xtrans);
    else
        PSHOW = find((1-Perm_Y2X_trans)<=PVALSHOW);
    end
    P3 = zeros(size(Perm_Y2X_trans));
    P3(PSHOW) = 1;
    Show_GC_res_subfun(Y2X_trans,P3,NROI,colormapshow,enhanind,outdir,'CaSCNresY2X_trans_perm',width1,width2);
    
    if PVALT==3
        PSHOW = find((1-Perm_FX2Y)<=PVALSHOW_PermFY2X);
    else
        PSHOW = find((1-Perm_FX2Y)<=PVALSHOW);
    end
    P3 = zeros(size(Perm_FX2Y));
    P3(PSHOW) = 1;
    Show_GC_res_subfun(FX2Y,P3,NROI,colormapshow,enhanind,outdir,'CaSCNresFY2X_perm',width1,width2);
end

end