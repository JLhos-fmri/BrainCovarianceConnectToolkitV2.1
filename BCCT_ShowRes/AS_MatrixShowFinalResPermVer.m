function AS_MatrixShowFinalResPermVer(enhanind)
if nargin<1
    enhanind = [];
end
[mpat,mnam,mext] = fileparts(which('AS_MatrixShowFinalResPermVer.m'));
load(fullfile(mpat,'gray.mat'))
colormapshow0 = AFNICOLORMAP(128);
% colormapshowO = AFNICOLORMAP(12);
colormapshow(1:64,:) = colormapshow0(1:64,:);
colormapshow(66:129,:) = colormapshow0(65:end,:);
colormapshow(65,:) = [1 1 1];

indir = uigetdir(pwd,'StatRes dir');
load(fullfile(indir,'PermPval.mat'));
load(fullfile(indir,'SetUpparameter.mat'));
outdir = uigetdir(pwd,'Outdir');

indir1 = Parameter.Inputdir1;
indir2 = Parameter.Inputdir2;
ShowRes1 = load(fullfile(indir1,'R_Pres.mat'));
ShowRes2 = load(fullfile(indir2,'R_Pres.mat'));


Hsize = get(0,'ScreenSize');
Bsize = min(Hsize(3),Hsize(4))*0.9;
Bsize2 = Bsize/3;
% POS1 = [100,100+Bsize2*2,Bsize2,Bsize2];
POS2 = [10,10+Bsize2*1,Bsize2*3,Bsize2];
% POS3 = [100,100+Bsize2*0,Bsize2,Bsize2];

ShowRes1.R(1:size(ShowRes1.R,1)+1:end) = 0;
ShowRes1.Z(1:size(ShowRes1.Z,1)+1:end) = 0;
ShowRes1.P(1:size(ShowRes1.P,1)+1:end) = 1;
ShowRes2.R(1:size(ShowRes2.R,1)+1:end) = 0;
ShowRes2.Z(1:size(ShowRes2.Z,1)+1:end) = 0;
ShowRes2.P(1:size(ShowRes2.P,1)+1:end) = 1;
NROI = size(ShowRes1.R,1);
enhanind = unique([0,NROI,enhanind]);

R = ShowRes1.R;
P = ShowRes1.P;
Z = ShowRes1.Z;
Pind005 = find(P<0.05);
PindBonf = find(P<0.05/((size(R,1)*(size(R,2)-1))/2));
P1 = zeros(size(P));
P1(Pind005) = 1;
P2 = zeros(size(P));
P2(PindBonf) = 1;
H = figure('pos',POS2,'name','Group 1: R');
ax = axes('parent',H,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
image(R,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end
% H = figure('pos',POS2,'name','Group 1: R(p<0.05)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
image(R.*P1,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end

for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end
% H = figure('pos',POS2,'name','Group 1: R(p<0.05,bonf corrected)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(R.*P2,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end

saveas(H,[outdir,filesep,'Group1R.fig'])
set(H,'PaperPositionMode','manual');
set(H,'PaperUnits','inch')
set(H,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H,[outdir,filesep,'Group1R.tif'],'-dtiff','-r300')
%%
H = figure('pos',POS2,'name','Group 1: Z');
ax = axes('parent',H,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
image(R,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end
% H = figure('pos',POS2,'name','Group 1: Z(p<0.05)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
image(R.*P1,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end

for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end
% H = figure('pos',POS2,'name','Group 1: Z(p<0.05,bonf corrected)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(R.*P2,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end

saveas(H,[outdir,filesep,'Group1Z.fig'])
set(H,'PaperPositionMode','manual');
set(H,'PaperUnits','inch')
set(H,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H,[outdir,filesep,'Group1Z.tif'],'-dtiff','-r300')
%%
R = ShowRes2.R;
P = ShowRes2.P;
Z = ShowRes2.Z;
Pind005 = find(P<0.05);
PindBonf = find(P<0.05/((size(R,1)*(size(R,2)-1))/2));
P1 = zeros(size(P));
P1(Pind005) = 1;
P2 = zeros(size(P));
P2(PindBonf) = 1;
H = figure('pos',POS2,'name','Group 2: R');
ax = axes('parent',H,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
image(R,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end
% H = figure('pos',POS2,'name','Group 2: R(p<0.05)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
image(R.*P1,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end

for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end
% H = figure('pos',POS2,'name','Group 2: R(p<0.05,bonf corrected)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(R.*P2,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end
saveas(H,[outdir,filesep,'Group2R.fig'])
set(H,'PaperPositionMode','manual');
set(H,'PaperUnits','inch')
set(H,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H,[outdir,filesep,'Group2R.tif'],'-dtiff','-r300')
%%
H = figure('pos',POS2,'name','Group 2: Z');
ax = axes('parent',H,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
image(R,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end
% H = figure('pos',POS2,'name','Group 2: Z(p<0.05)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
image(R.*P1,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end

for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end
% H = figure('pos',POS2,'name','Group 2: Z(p<0.05,bonf corrected)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(R.*P2,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end

saveas(H,[outdir,filesep,'Group2Z.fig'])

set(H,'PaperPositionMode','manual');
set(H,'PaperUnits','inch')
set(H,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H,[outdir,filesep,'Group2Z.tif'],'-dtiff','-r300')

THRBONF = 0.05/((NROI*(NROI-1))/2);
%%
PvalShow = Pval;
ZPvalShow = ZPval;
PvalShow(1:NROI+1:end) = 0.5;
ZPvalShow(1:NROI+1:end) = 0.5;

PvalShow(Pval<0.5) = PvalShow(Pval<0.5)-1;
ZPvalShow(ZPval<0.5) = ZPvalShow(ZPval<0.5)-1;
PvalShow(1:NROI+1:end) = 0;
ZPvalShow(1:NROI+1:end) = 0;
H = figure('pos',POS2,'name','Statistical P(R)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
image(PvalShow,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end


H2 = figure('pos',POS2,'name','Statistical P(Z)');
ax = axes('parent',H2,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
image(ZPvalShow,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end


%
k = 32;
mk = 16;
clear colormapshowP1new
colormapshowP1 = AFNICOLORMAP(k);
% indnum = floor((0.95/0.05)*6);
colormapshowP1new(1:mk,:) = colormapshowP1(1:mk,:);
colormapshowP1new(1+mk:mk*3+1,:) = 1;
colormapshowP1new(mk*3+1+1:mk*4+1,:) = colormapshowP1(mk+1:k,:);

PvalShow2 = PvalShow;
PvalShow2(PvalShow<0.95&PvalShow>-0.95) = 0;
PvalShow2(PvalShow>=0.95) = PvalShow2(PvalShow>=0.95)-0.9;
PvalShow2(PvalShow<=-0.95) = PvalShow2(PvalShow<=-0.95)+0.9;
ZPvalShow2 = ZPvalShow;
ZPvalShow2(ZPvalShow<0.95&ZPvalShow>-0.95) = 0;
ZPvalShow2(ZPvalShow>=0.95) = ZPvalShow2(ZPvalShow>=0.95)-0.9;
ZPvalShow2(ZPvalShow<=-0.95) = ZPvalShow2(ZPvalShow<=-0.95)+0.9;
% H = figure('pos',POS2,'name','Statistical P<0.05(R)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
image(PvalShow2,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshowP1new);
set(ax,'Clim',[-0.1,0.1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end

%%
% H = figure('pos',POS2,'name','Statistical P<0.05(Z)');
ax = axes('parent',H2,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
image(ZPvalShow2,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshowP1new);
set(ax,'Clim',[-0.1,0.1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end

%
PvalShow3 = PvalShow;
THRV2 = 1-THRBONF;
THRV3 = 1-THRBONF*2;
PvalShow3(PvalShow<1-THRBONF&PvalShow>THRBONF-1) = 0;
PvalShow2(PvalShow>=THRV2) = PvalShow2(PvalShow>=THRV2)-THRV3;
PvalShow2(PvalShow<=-THRV2) = PvalShow2(PvalShow<=-THRV2)+THRV3;

ZPvalShow3 = ZPvalShow;
ZPvalShow3(ZPvalShow<1-THRBONF&ZPvalShow>THRBONF-1) = 0;
ZPvalShow2(ZPvalShow>=THRV2) = ZPvalShow2(ZPvalShow>=THRV2)-THRV3;
ZPvalShow2(ZPvalShow<=-THRV2) = ZPvalShow2(ZPvalShow<=-THRV2)+THRV3;

% H = figure('pos',POS2,'name','Statistical P<0.05 Bonf corrected(R)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(PvalShow3,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshowP1new);
set(ax,'Clim',[-1,1]*THRBONF*2)
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end

% H = figure('pos',POS2,'name','Statistical P<0.05 Bonf corrected(z)');
ax = axes('parent',H2,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(ZPvalShow3,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshowP1new);
set(ax,'Clim',[-1,1]*THRBONF*2)
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',0.5,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',0.5,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',3,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',3,'parent',ax);
end

saveas(H,[outdir,filesep,'Stat_R.fig'])

set(H,'PaperPositionMode','manual');
set(H,'PaperUnits','inch')
set(H,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H,[outdir,filesep,'Stat_R.tif'],'-dtiff','-r300')
saveas(H2,[outdir,filesep,'Stat_Z.fig'])
set(H2,'PaperPositionMode','manual');
set(H2,'PaperUnits','inch')
set(H2,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H2,[outdir,filesep,'Stat_Z.tif'],'-dtiff','-r300')
end