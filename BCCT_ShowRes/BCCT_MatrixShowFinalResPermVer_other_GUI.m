function BCCT_MatrixShowFinalResPermVer_other_GUI(Pval,outdir,enhanind,Coltype,PVALT,Pvalused,orders,outname)
[mpat,mnam,mext] = fileparts(which('BCCT_MatrixShowFinalResPermVer_other_GUI.m'));
load(fullfile(mpat,'MatColormap','mycolormap.mat'));
load(fullfile(mpat,'gray.mat'))
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

Hsize = get(0,'ScreenSize');
Bsize = min(Hsize(3),Hsize(4))*0.9;
Bsize2 = Bsize;
POS3 = [10,10,Bsize2,Bsize2];

Pval = Pval(orders,orders);
NROI = size(Pval,1);
if NROI>50
    width1 = 0.01;
    width2 = 1;
else
    width1 = 0.5;
    width2 = 3;
end
enhanind = unique([0,NROI,enhanind]);

%%
PvalShow = Pval;
% ZPvalShow = ZPval;
PvalShow(1:NROI+1:end) = 0.5;
% ZPvalShow(1:NROI+1:end) = 0.5;

%%
H_r = figure('pos',POS3,'name',outname);
ax_r = axes('parent',H_r,'unit','norm','pos',[0.1,0.1,0.8,0.8]);
% H_z = figure('pos',POS3,'name','Statistical P(R) Setup');
% ax_z = axes('parent',H_z,'unit','norm','pos',[0.1,0.1,0.8,0.8]);
Nshowtempnud = ones(NROI);
Nshowtempnud(1:NROI+1:end) = 0;
if PVALT==1
    PVALSHOW = Pvalused;
%     PVALSHOWZ = PVALSHOW;
elseif PVALT==2
    PVALSHOW = 1/(NROI*(NROI-1)/2);
%     PVALSHOWZ = PVALSHOW;
elseif PVALT==3
    Pusedinfo = Pval(find(Nshowtempnud));
%     PusedinfoZ = ZPval(find(Nshowtempnud));
    [h pi] = fdr(Pusedinfo, Pvalused);
%     [hz piz] = fdr(PusedinfoZ, Pvalused);
    hind = find(h);
    hindz = find(h);
    if ~isempty(hind)
        PVALSHOW = max(pi(hind));
    else
        PVALSHOW = 1e-30;
    end
else
    PVALSHOW = Pvalused/(NROI*(NROI-1)/2);
end
if PVALSHOW<1e-6
    PVALSHOW1 = 1e-6;
else
    PVALSHOW1 = PVALSHOW;
end

PvalShowTEMP = Pval;
PvalShowTEMP(1:NROI+1:end) = 0.5;
% ZPvalShowTEMP(1:NROI+1:end) = 0.5;
A2 = zeros(size(PvalShowTEMP));
A2(PvalShowTEMP<=PVALSHOW1) = PvalShowTEMP(PvalShowTEMP<=PVALSHOW1)-1;
A2(PvalShowTEMP>=1-PVALSHOW1) = PvalShowTEMP(PvalShowTEMP>=1-PVALSHOW1);

image(A2,'parent',ax_r,'CDataMapping','scaled');
axis(ax_r,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax_r,'off');
colormap(ax_r,colormapshow);
set(ax_r,'Clim',[-1,1])
hold(ax_r,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax_r);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax_r);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax_r);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax_r);
end


saveas(H_r,[outdir,filesep,outname,'.fig'])
set(H_r,'PaperPositionMode','manual');
set(H_r,'PaperUnits','inch')
set(H_r,'Paperposition',[1,1,POS3(3)*3/300,POS3(4)*3/300]);
print(H_r,[outdir,filesep,outname,'.tif'],'-dtiff','-r300')
end