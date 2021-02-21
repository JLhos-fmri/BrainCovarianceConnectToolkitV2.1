function AS_MatrixFCshow(enhanind)
if nargin<1
    enhanind = [];
end
indir = uigetdir(pwd,'MatrixFC dir');
load(fullfile(indir,'R_Pres.mat'));
[mpat,mnam,mext] = fileparts(which('AS_MatrixFCshow.m'));
load(fullfile(mpat,'gray.mat'))
outdir = uigetdir(pwd,'Outdir');

R(1:size(R,1)+1:end) = 0;
Z(1:size(Z,1)+1:end) = 0;
P(1:size(P,1)+1:end) = 1;
NROI = size(R,1);
Hsize = get(0,'ScreenSize');
Bsize = min(Hsize(3),Hsize(4))*0.9;
Bsize2 = Bsize/3;
POS1 = [10,10+Bsize2*2,Bsize2*3,Bsize2];
POS2 = [10,10+Bsize2*1,Bsize2*3,Bsize2];
POS3 = [10,10+Bsize2*0,Bsize2*3,Bsize2];
H1 = figure('pos',POS1);
H2 = figure('pos',POS2);
H3 = figure('pos',POS3);
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

% colormapshow = AFNICOLORMAP(128);
% % colormapshowO = AFNICOLORMAP(12);
% colormapshow(64,:) = [1 1 1];
% colormapshow(65,:) = [1 1 1];

colormapshow0 = AFNICOLORMAP(128);
% colormapshowO = AFNICOLORMAP(12);
colormapshow(1:64,:) = colormapshow0(1:64,:);
colormapshow(66:129,:) = colormapshow0(65:end,:);
colormapshow(65,:) = [1 1 1];
enhanind = unique([0,NROI,enhanind]);

image(R,'parent',axH1_1,'CDataMapping','scaled');
axis(axH1_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH1_1,'off');
colormap(axH1_1,colormapshow);
set(axH1_1,'Clim',[-1,1])
hold(axH1_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',axH1_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',axH1_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',axH1_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',axH1_1);
end

image(R.*P1,'parent',axH1_2,'CDataMapping','scaled');
axis(axH1_2,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH1_2,'off');
colormap(axH1_2,colormapshow);
set(axH1_2,'Clim',[-1,1])
hold(axH1_2,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',axH1_2);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',axH1_2);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',axH1_2);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',axH1_2);
end

image(R.*P2,'parent',axH1_3,'CDataMapping','scaled');
axis(axH1_3,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH1_3,'off');
colormap(axH1_3,colormapshow);
set(axH1_3,'Clim',[-1,1])
hold(axH1_3,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',axH1_3);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',axH1_3);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',axH1_3);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',axH1_3);
end
%
image(Z,'parent',axH2_1,'CDataMapping','scaled');
axis(axH2_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH2_1,'off');
colormap(axH2_1,colormapshow);
set(axH2_1,'Clim',[-1,1]*maxZ)
hold(axH2_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',axH2_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',axH2_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',axH2_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',axH2_1);
end

image(Z.*P1,'parent',axH2_2,'CDataMapping','scaled');
axis(axH2_2,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH2_2,'off');
colormap(axH2_2,colormapshow);
set(axH2_2,'Clim',[-1,1]*maxZ)
hold(axH2_2,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',axH2_2);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',axH2_2);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',axH2_2);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',axH2_2);
end


image(Z.*P2,'parent',axH2_3,'CDataMapping','scaled');
axis(axH2_3,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH2_3,'off');
colormap(axH2_3,colormapshow);
set(axH2_3,'Clim',[-1,1]*maxZ)
hold(axH2_3,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',axH2_3);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',axH2_3);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',axH2_3);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',axH2_3);
end


%
% Colormgray = colormapshow0;
% Colormgray(1:64,:) = 1;
Colormgrayp = colormapshow0;
Colormgrayp(1:64,:) = 0.8;
% Colormgray(1,:) = [1 1 1];
% overthr = 0.1;
overthr = 0;
Pnew = 1-P;
Pnew(1:NROI+1:end) = 0;
image(Pnew,'parent',axH3_1,'CDataMapping','scaled');
axis(axH3_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH3_1,'off');
colormap(axH3_1,Colormgrayp);
set(axH3_1,'Clim',[0,1+overthr*0.5])
hold(axH3_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',axH3_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',axH3_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',axH3_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',axH3_1);
end

image(Pnew.*P1,'parent',axH3_2,'CDataMapping','scaled');
axis(axH3_2,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH3_2,'off');
colormap(axH3_2,Colormgrayp);
set(axH3_2,'Clim',[0.90,1+overthr*0.05])
hold(axH3_2,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',axH3_2);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',axH3_2);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',axH3_2);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',axH3_2);
end

image(Pnew.*P2,'parent',axH3_3,'CDataMapping','scaled');
axis(axH3_3,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axH3_3,'off');
colormap(axH3_3,Colormgrayp);
pbonfthr = 0.05/((size(R,1)*(size(R,2)-1))/2);
set(axH3_3,'Clim',[1-pbonfthr*2,1+overthr*(pbonfthr)]);
hold(axH3_3,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',axH3_3);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',axH3_3);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',axH3_3);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',axH3_3);
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
end