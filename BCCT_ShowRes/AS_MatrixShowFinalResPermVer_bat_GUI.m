function AS_MatrixShowFinalResPermVer_bat_GUI(indir,outdir,enhanind,Coltype,PVALT,Pvalused,orders)

[mpat,mnam,mext] = fileparts(which('AS_MatrixShowFinalResPermVer_bat_GUI.m'));
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

% indir = uigetdir(pwd,'StatRes dir');
% load(fullfile(indir,'PermPval.mat'));
load(indir)
% load(fullfile(indir,'SetUpparameter.mat'));
% outdir = uigetdir(pwd,'Outdir');

% indir1 = Parameter.Inputdir1;
% indir2 = Parameter.Inputdir2;
% ShowRes1 = load(fullfile(indir1,'R_Pres.mat'));
% ShowRes2 = load(fullfile(indir2,'R_Pres.mat'));

Hsize = get(0,'ScreenSize');
Bsize = min(Hsize(3),Hsize(4))*0.9;
Bsize2 = Bsize/3;
POS2 = [10,10+Bsize2*1,Bsize2*3,Bsize2];
POS3 = [10,10+Bsize2*1,Bsize2*1,Bsize2];

% ShowRes1.R(1:size(ShowRes1.R,1)+1:end) = 0;
% ShowRes1.Z(1:size(ShowRes1.Z,1)+1:end) = 0;
% ShowRes1.P(1:size(ShowRes1.P,1)+1:end) = 1;
% ShowRes2.R(1:size(ShowRes2.R,1)+1:end) = 0;
% ShowRes2.Z(1:size(ShowRes2.Z,1)+1:end) = 0;
% ShowRes2.P(1:size(ShowRes2.P,1)+1:end) = 1;
Pval = P_mat;
ZPval = P_mat;
Pval = Pval(orders,orders);
ZPval = ZPval(orders,orders);
NROI = size(Pval,1);
if NROI>50
    width1 = 0.01;
    width2 = 1;
else
    width1 = 0.5;
    width2 = 3;
end
enhanind = unique([0,NROI,enhanind]);

%% Result of Group1 R
% R = ShowRes1.R;
% P = ShowRes1.P;
% Z = ShowRes1.Z;
% Pind005 = find(P<0.05);
% PindBonf = find(P<0.05/((size(R,1)*(size(R,2)-1))/2));
% P1 = zeros(size(P));
% P1(Pind005) = 1;
% P2 = zeros(size(P));
% P2(PindBonf) = 1;
% H = figure('pos',POS2,'name','Group 1: R');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
% image(R,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% % H = figure('pos',POS2,'name','Group 1: R(p<0.05)');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
% image(R.*P1,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% 
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% % H = figure('pos',POS2,'name','Group 1: R(p<0.05,bonf corrected)');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
% image(R.*P2,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% 
% saveas(H,[outdir,filesep,'Group1R.fig'])
% set(H,'PaperPositionMode','manual');
% set(H,'PaperUnits','inch')
% set(H,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
% print(H,[outdir,filesep,'Group1R.tif'],'-dtiff','-r300')

%% Result of Group1 Z
% H = figure('pos',POS2,'name','Group 1: Z');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
% image(R,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% % H = figure('pos',POS2,'name','Group 1: Z(p<0.05)');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
% image(R.*P1,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% 
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% % H = figure('pos',POS2,'name','Group 1: Z(p<0.05,bonf corrected)');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
% image(R.*P2,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% 
% saveas(H,[outdir,filesep,'Group1Z.fig'])
% set(H,'PaperPositionMode','manual');
% set(H,'PaperUnits','inch')
% set(H,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
% print(H,[outdir,filesep,'Group1Z.tif'],'-dtiff','-r300')
%% Result of Group 2
% R = ShowRes2.R;
% P = ShowRes2.P;
% Z = ShowRes2.Z;
% Pind005 = find(P<0.05);
% PindBonf = find(P<0.05/((size(R,1)*(size(R,2)-1))/2));
% P1 = zeros(size(P));
% P1(Pind005) = 1;
% P2 = zeros(size(P));
% P2(PindBonf) = 1;
% H = figure('pos',POS2,'name','Group 2: R');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
% image(R,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% % H = figure('pos',POS2,'name','Group 2: R(p<0.05)');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
% image(R.*P1,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% 
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% % H = figure('pos',POS2,'name','Group 2: R(p<0.05,bonf corrected)');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
% image(R.*P2,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% saveas(H,[outdir,filesep,'Group2R.fig'])
% set(H,'PaperPositionMode','manual');
% set(H,'PaperUnits','inch')
% set(H,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
% print(H,[outdir,filesep,'Group2R.tif'],'-dtiff','-r300')
%% Result of Group2 Z
% H = figure('pos',POS2,'name','Group 2: Z');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
% image(R,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% % H = figure('pos',POS2,'name','Group 2: Z(p<0.05)');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3+1/3,0.1,0.8/3,0.8]);
% image(R.*P1,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% 
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% % H = figure('pos',POS2,'name','Group 2: Z(p<0.05,bonf corrected)');
% ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
% image(R.*P2,'parent',ax,'CDataMapping','scaled');
% axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
% axis(ax,'off');
% colormap(ax,colormapshow);
% set(ax,'Clim',[-1,1])
% hold(ax,'on');
% for i = 1:NROI
%     plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
%     plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
% end
% for i = 1:length(enhanind)
%     plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
%     plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
% end
% 
% saveas(H,[outdir,filesep,'Group2Z.fig'])
% 
% set(H,'PaperPositionMode','manual');
% set(H,'PaperUnits','inch')
% set(H,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
% print(H,[outdir,filesep,'Group2Z.tif'],'-dtiff','-r300')


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
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
end

%%
H2 = figure('pos',POS2,'name','Statistical P(Z)');
ax = axes('parent',H2,'unit','norm','pos',[0.1/3,0.1,0.8/3,0.8]);
image(ZPvalShow,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshow);
set(ax,'Clim',[-1,1])
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
end

%%
%
k = 32;
mk = 16;
clear colormapshowP1new
colormapshowP1 = AFNICOLORMAP(k);
% indnum = floor((0.95/0.05)*6);
colormapshowP1new(1:mk,:) = colormapshowP1(1:mk,:);
colormapshowP1new(1+mk:mk*3+1,:) = 0.8;
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
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
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
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
end

%
PvalShow3 = PvalShow;
THRV2 = 1-THRBONF;
THRV3 = 1-THRBONF*2;
PvalShow3(PvalShow<1-THRBONF&PvalShow>THRBONF-1) = 0;
PvalShow3(PvalShow>=THRV2) = PvalShow3(PvalShow>=THRV2)-THRV3;
PvalShow3(PvalShow<=-THRV2) = PvalShow3(PvalShow<=-THRV2)+THRV3;

ZPvalShow3 = ZPvalShow;
ZPvalShow3(ZPvalShow<1-THRBONF&ZPvalShow>THRBONF-1) = 0;
ZPvalShow3(ZPvalShow>=THRV2) = ZPvalShow3(ZPvalShow>=THRV2)-THRV3;
ZPvalShow3(ZPvalShow<=-THRV2) = ZPvalShow3(ZPvalShow<=-THRV2)+THRV3;

% H = figure('pos',POS2,'name','Statistical P<0.05 Bonf corrected(R)');
ax = axes('parent',H,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(PvalShow3,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshowP1new);
set(ax,'Clim',[-1,1]*THRBONF*2)
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
end

ax = axes('parent',H2,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(ZPvalShow3,'parent',ax,'CDataMapping','scaled');
axis(ax,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax,'off');
colormap(ax,colormapshowP1new);
set(ax,'Clim',[-1,1]*THRBONF*2)
hold(ax,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax);
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

%%
H_r = figure('pos',POS3,'name','Statistical P(R) Setup');
ax_r = axes('parent',H_r,'unit','norm','pos',[0.1,0.1,0.8,0.8]);
H_z = figure('pos',POS3,'name','Statistical P(R) Setup');
ax_z = axes('parent',H_z,'unit','norm','pos',[0.1,0.1,0.8,0.8]);
Nshowtempnud = tril(ones(NROI),-1);
if PVALT==1
    PVALSHOW = Pvalused;
    PVALSHOWZ = PVALSHOW;
elseif PVALT==2
    PVALSHOW = 1/(NROI*(NROI-1)/2);
    PVALSHOWZ = PVALSHOW;
elseif PVALT==3
    Pusedinfo = Pval(find(Nshowtempnud));
    PusedinfoZ = ZPval(find(Nshowtempnud));
    [h pi] = fdr(Pusedinfo, Pvalused);
    [hz piz] = fdr(PusedinfoZ, Pvalused);
    hind = find(h);
    hindz = find(h);
    if ~isempty(hind)
        PVALSHOW = max(pi(hind));
    else
        PVALSHOW = 1e-30;
    end
    if ~isempty(hindz)
        PVALSHOWZ = max(piz(hindz));
    else
        PVALSHOWZ = 1e-30;
    end
else
    PVALSHOW = Pvalused/(NROI*(NROI-1)/2);
    PVALSHOWZ = PVALSHOW;
end
if PVALSHOW<1e-6
    PVALSHOW1 = 1e-6;
else
    PVALSHOW1 = PVALSHOW;
end

if PVALSHOWZ<1e-6
    PVALSHOW2 = 1e-6;
else
    PVALSHOW2 = PVALSHOWZ;
end

PvalShowTEMP = Pval;
ZPvalShowTEMP = ZPval;
PvalShowTEMP(1:NROI+1:end) = 0.5;
ZPvalShowTEMP(1:NROI+1:end) = 0.5;
PvalShowTEMP(Pval>0.5) = 1-Pval((Pval>0.5));
ZPvalShowTEMP(ZPval>0.5) = 1-ZPval((ZPval>0.5));

PSHOW = find(PvalShowTEMP<=PVALSHOW1);
P3 = zeros(size(Pval));
P3(PSHOW) = 1;

PSHOWZ = find(ZPvalShowTEMP<=PVALSHOW2);
P4 = zeros(size(Pval));
P4(PSHOWZ) = 1;
PvalShowS = PvalShow;
ZPvalShowS = ZPvalShow;

THRV2 = 1-PVALSHOW1;
THRV3 = 1-PVALSHOW1*2;

THRV2Z = 1-PVALSHOW2;
THRV3Z = 1-PVALSHOW2*2;

PvalShowS(PvalShow<1-PVALSHOW&PvalShow>PVALSHOW-1) = 0;
PvalShowS(PvalShow>=THRV2) = PvalShowS(PvalShow>=THRV2)-THRV3;
PvalShowS(PvalShow<=-THRV2) = PvalShowS(PvalShow<=-THRV2)+THRV3;

ZPvalShowS(ZPvalShow<1-PVALSHOWZ&ZPvalShow>PVALSHOWZ-1) = 0;
ZPvalShowS(ZPvalShow>=THRV2Z) = ZPvalShowS(ZPvalShow>=THRV2Z)-THRV3Z;
ZPvalShowS(ZPvalShow<=-THRV2Z) = ZPvalShowS(ZPvalShow<=-THRV2Z)+THRV3Z;
%
% ax = axes('parent',H_r,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(PvalShowS,'parent',ax_r,'CDataMapping','scaled');
axis(ax_r,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax_r,'off');
colormap(ax_r,colormapshowP1new);
set(ax_r,'Clim',[-1,1]*THRBONF*2)
hold(ax_r,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax_r);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax_r);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax_r);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax_r);
end

% ax = axes('parent',H2,'unit','norm','pos',[0.1/3+2/3,0.1,0.8/3,0.8]);
image(ZPvalShowS,'parent',ax_z,'CDataMapping','scaled');
axis(ax_z,[-0.5,NROI+1.5,-0.5,NROI+1.5]);
axis(ax_z,'off');
colormap(ax_z,colormapshowP1new);
set(ax_z,'Clim',[-1,1]*THRBONF*2)
hold(ax_z,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',ax_z);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',ax_z);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',ax_z);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',ax_z);
end
saveas(H_r,[outdir,filesep,'Stat_R_setup.fig'])
set(H_r,'PaperPositionMode','manual');
set(H_r,'PaperUnits','inch')
set(H_r,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H_r,[outdir,filesep,'Stat_R_setup.tif'],'-dtiff','-r300')
saveas(H_z,[outdir,filesep,'Stat_Z_setup.fig'])
set(H_z,'PaperPositionMode','manual');
set(H_z,'PaperUnits','inch')
set(H_z,'Paperposition',[1,1,POS2(3)*3/300,POS2(4)*3/300]);
print(H_z,[outdir,filesep,'Stat_Z_setup.tif'],'-dtiff','-r300')
end