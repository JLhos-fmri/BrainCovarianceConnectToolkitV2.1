function AS_MatrixFCshow_other_GUI(R,outdir,enhanind,Coltype,PVALT,Pvalused,orders,dof)

% load(indir);
[mpat,mnam,mext] = fileparts(which('AS_MatrixFCshow_other_GUI.m'));
load(fullfile(mpat,'gray.mat'));
load(fullfile(mpat,'MatColormap','mycolormap.mat'));

R(1:size(R,1)+1:end) = 0;
R = R(orders,orders);
dof2 = [];
[Z, P] = AS_TFRtoZ(R,'R',dof,dof2);
Z(1:size(Z,1)+1:end) = 0;
P(1:size(P,1)+1:end) = 1;
Z(isnan(Z)) = 0;
P(isnan(P)) = 1;
Z(isinf(Z)) = 0;
P(isinf(P)) = 1;
Nshowtempnud = tril(ones(size(P)),-1);

NROI = size(R,1);
if NROI>50
    width1 = 0.01;
    width2 = 1;
else
    width1 = 0.5;
    width2 = 3;
end
%%
if PVALT==1
    PVALSHOW = Pvalused;
elseif PVALT==2
    PVALSHOW = 1/(NROI*(NROI-1)/2);
elseif PVALT==3
    Pusedinfo = P(find(Nshowtempnud));
    [h pi] = fdr(Pusedinfo, Pvalused);
    hind = find(h);
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
PSHOW = find(P<=PVALSHOW);
P3 = zeros(size(P));
P3(PSHOW) = 1;

Hsize = get(0,'ScreenSize');
Bsize = min(Hsize(3),Hsize(4))*0.9;
Bsize2 = Bsize/3;
POS1 = [10,10+Bsize2*2,Bsize2*3,Bsize2];
POS2 = [10,10+Bsize2*1,Bsize2*3,Bsize2];
POS3 = [10,10+Bsize2*0,Bsize2*3,Bsize2];
POSm1 = [20,20+Bsize2*2,Bsize2,Bsize2];
POSm2 = [20,20+Bsize2*1,Bsize2,Bsize2];
POSm3 = [20,20+Bsize2*0,Bsize2,Bsize2];
H1 = figure('pos',POS1);
H2 = figure('pos',POS2);
H3 = figure('pos',POS3);
Hm1 = figure('pos',POSm1);
Hm2 = figure('pos',POSm2);
Hm3 = figure('pos',POSm3);
maxZ = max(abs(Z(:)));
Pind005 = find(P<0.05);
PindBonf = find(P<0.05/((size(R,1)*(size(R,2)-1))/2));
P1 = zeros(size(P));
P1(Pind005) = 1;
P2 = zeros(size(P));
P2(PindBonf) = 1;
axH1_1 = axes('parent',H1,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8],'CLim',[-1 1]);
axH1_2 = axes('parent',H1,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8],'CLim',[-1 1]);
axH1_3 = axes('parent',H1,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8],'CLim',[-1 1]);
axH2_1 = axes('parent',H2,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8],'CLim',[-1 1]*maxZ);
axH2_2 = axes('parent',H2,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8],'CLim',[-1 1]*maxZ);
axH2_3 = axes('parent',H2,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8],'CLim',[-1 1]*maxZ);
axH3_1 = axes('parent',H3,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8],'CLim',[0 1]);
axH3_2 = axes('parent',H3,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8],'CLim',[0 1]);
axH3_3 = axes('parent',H3,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8],'CLim',[0 1]);

axHm1_1 = axes('parent',Hm1,'unit','norm','pos',[0.1 0.1 0.8 0.8],'CLim',[-1 1]);
axHm2_1 = axes('parent',Hm2,'unit','norm','pos',[0.1 0.1 0.8 0.8],'CLim',[-1 1]*maxZ);
axHm3_1 = axes('parent',Hm3,'unit','norm','pos',[0.1 0.1 0.8 0.8],'CLim',[0 1]);
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

image(R,'parent',axH1_1,'CDataMapping','scaled');
axis(axH1_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH1_1,'off');
colormap(axH1_1,colormapshow);
set(axH1_1,'Clim',[-1,1])
hold(axH1_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axH1_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axH1_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axH1_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axH1_1);
end

image(R.*P1,'parent',axH1_2,'CDataMapping','scaled');
axis(axH1_2,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH1_2,'off');
colormap(axH1_2,colormapshow);
set(axH1_2,'Clim',[-1,1])
hold(axH1_2,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axH1_2);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axH1_2);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axH1_2);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axH1_2);
end

image(R.*P2,'parent',axH1_3,'CDataMapping','scaled');
axis(axH1_3,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH1_3,'off');
colormap(axH1_3,colormapshow);
set(axH1_3,'Clim',[-1,1])
hold(axH1_3,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axH1_3);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axH1_3);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axH1_3);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axH1_3);
end
%
image(Z,'parent',axH2_1,'CDataMapping','scaled');
axis(axH2_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH2_1,'off');
colormap(axH2_1,colormapshow);
set(axH2_1,'Clim',[-1,1]*maxZ)
hold(axH2_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axH2_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axH2_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axH2_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axH2_1);
end

image(Z.*P1,'parent',axH2_2,'CDataMapping','scaled');
axis(axH2_2,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH2_2,'off');
colormap(axH2_2,colormapshow);
set(axH2_2,'Clim',[-1,1]*maxZ)
hold(axH2_2,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axH2_2);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axH2_2);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axH2_2);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axH2_2);
end


image(Z.*P2,'parent',axH2_3,'CDataMapping','scaled');
axis(axH2_3,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH2_3,'off');
colormap(axH2_3,colormapshow);
set(axH2_3,'Clim',[-1,1]*maxZ)
hold(axH2_3,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axH2_3);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axH2_3);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axH2_3);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axH2_3);
end


%
Colormgray = colormapshow0;
Colormgray(1:64,:) = 0.8;
% Colormgray(1,:) = [1 1 1];
% overthr = 0.1;
overthr = 0;
Pnew = 1-P;
Pnew(1:NROI+1:end) = 0;
image(Pnew,'parent',axH3_1,'CDataMapping','scaled');
axis(axH3_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH3_1,'off');
colormap(axH3_1,Colormgray);
set(axH3_1,'Clim',[0,1+overthr*0.5])
hold(axH3_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axH3_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axH3_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axH3_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axH3_1);
end

image(Pnew.*P1,'parent',axH3_2,'CDataMapping','scaled');
axis(axH3_2,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH3_2,'off');
colormap(axH3_2,Colormgray);
set(axH3_2,'Clim',[0.90,1+overthr*0.05])
hold(axH3_2,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axH3_2);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axH3_2);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axH3_2);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axH3_2);
end

image(Pnew.*P2,'parent',axH3_3,'CDataMapping','scaled');
axis(axH3_3,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH3_3,'off');
colormap(axH3_3,Colormgray);
pbonfthr = 0.05/((size(R,1)*(size(R,2)-1))/2);
set(axH3_3,'Clim',[1-pbonfthr*2,1+overthr*(pbonfthr)]);
hold(axH3_3,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axH3_3);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axH3_3);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axH3_3);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axH3_3);
end

% saveas(h,[OutputdirFig,filesep,savename,'.fig']);
saveas(H1,[outdir,filesep,'R.fig'])
saveas(H2,[outdir,filesep,'Z.fig'])
saveas(H3,[outdir,filesep,'p.fig'])
% fpath=fullfile(pathname,filename);
% [pathstr, name, ext] = fileparts(fpath);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inch');
% set(gcf,'Paperposition',[1 1 EC.img.width/EC.img.dpi EC.img.height/EC.img.dpi]);
set(H1,'PaperPositionMode','manual');
set(H1,'PaperUnits','inch')
set(H1,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H1,[outdir,filesep,'R.tif'],'-dtiff','-r300')

set(H2,'PaperPositionMode','manual');
set(H2,'PaperUnits','inch')
set(H2,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H2,[outdir,filesep,'Z.tif'],'-dtiff','-r300')

set(H3,'PaperPositionMode','manual');
set(H3,'PaperUnits','inch')
set(H3,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H3,[outdir,filesep,'P.tif'],'-dtiff','-r300')

%%


image(R.*P3,'parent',axHm1_1,'CDataMapping','scaled');
axis(axHm1_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axHm1_1,'off');
colormap(axHm1_1,colormapshow);
set(axHm1_1,'Clim',[-1,1])
hold(axHm1_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axHm1_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axHm1_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axHm1_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axHm1_1);
end

image(Z.*P3,'parent',axHm2_1,'CDataMapping','scaled');
axis(axHm2_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axHm2_1,'off');
colormap(axHm2_1,colormapshow);
set(axHm2_1,'Clim',[-1,1]*maxZ)
hold(axHm2_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axHm2_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axHm2_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axHm2_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axHm2_1);
end

image(Pnew.*P3,'parent',axHm3_1,'CDataMapping','scaled');
axis(axHm3_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axHm3_1,'off');
colormap(axHm3_1,Colormgray);
% PSHOW = 0.05/((size(R,1)*(size(R,2)-1))/2);
set(axHm3_1,'Clim',[1-PVALSHOW1*2,1+overthr*(PVALSHOW1)]);
hold(axHm3_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axHm3_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axHm3_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axHm3_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axHm3_1);
end

saveas(Hm1,[outdir,filesep,'R_setup.fig'])
saveas(Hm2,[outdir,filesep,'Z_setup.fig'])
saveas(Hm3,[outdir,filesep,'p_setup.fig'])
% fpath=fullfile(pathname,filename);
% [pathstr, name, ext] = fileparts(fpath);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inch');
% set(gcf,'Paperposition',[1 1 EC.img.width/EC.img.dpi EC.img.height/EC.img.dpi]);
set(Hm1,'PaperPositionMode','manual');
set(Hm1,'PaperUnits','inch')
set(Hm1,'Paperposition',[1,1,POSm1(3)*3/300,POSm1(4)*3/300]);
print(Hm1,[outdir,filesep,'R_setup.tif'],'-dtiff','-r300')

set(Hm2,'PaperPositionMode','manual');
set(Hm2,'PaperUnits','inch')
set(Hm2,'Paperposition',[1,1,POSm2(3)*3/300,POSm2(4)*3/300]);
print(Hm2,[outdir,filesep,'Z_setup.tif'],'-dtiff','-r300')

set(Hm3,'PaperPositionMode','manual');
set(Hm3,'PaperUnits','inch')
set(Hm3,'Paperposition',[1,1,POSm3(3)*3/300,POSm3(4)*3/300]);
print(Hm3,[outdir,filesep,'P_setup.tif'],'-dtiff','-r300')
end