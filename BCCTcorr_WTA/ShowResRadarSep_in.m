function ShowResRadarSep_in(Outputdir)
% Outputdir = uigetdir(pwd,'Output Directory Selection');
OutputdirFig = [Outputdir,filesep,'ResRadarFigure'];
if isempty(dir(OutputdirFig))
    mkdir(OutputdirFig)
end
% Outputdir = 'D:\testdata\AS_sub2cor1\testdata\';
LabedValmat = fullfile(Outputdir,'LabedVal.mat');
RealComputemat = fullfile(Outputdir,'RealCompute.mat');
load(RealComputemat)
computevalmat = fullfile(Outputdir,'computeval.mat');
load(LabedValmat)
load(computevalmat)
[IX IY] = sort(maxind);
for i = 1:seednum
    SX(i,1) = length(find(IX==i));
    SXpercent(i,1) = SX(i,1)/length(IX);
end
maxpercent = max(SXpercent);

pathmfile = which('AS_WTA_GUI.m');
[pat nam ext] = fileparts(pathmfile);
load(fullfile(pat,'SEEDCOLOR.mat'))
load(fullfile(pat,'JETcode.mat'));
Hsize = get(0,'ScreenSize');
Bsize = min(Hsize(3),Hsize(4));
BsizeUsed = floor(Bsize*0.6);
POS = [ceil((Hsize(3)-BsizeUsed)/2),ceil((Hsize(4)-BsizeUsed)/2),BsizeUsed,BsizeUsed];

h = figure('pos',POS);

ax = axes('parent',h,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax,'on');
% for i = 1:5
%     CirclePlot(i/5,ax);
% end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)],[0,YShow(i,1)],'k--','parent',ax);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
    for j = 1:5
        plot([XShow2(i,1)*j/5,XShow2(i+1,1)*j/5],...
            [YShow2(i,1)*j/5,YShow2(i+1,1)*j/5],...
            '--','parent',ax,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)],...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)],...
        'parent',ax,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i),YShow2(i,1)*SXpercent2(i),'parent',ax,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax,[-1,1,-1,1])
title('Percent for Each Seed ROI')
saveas(h,[OutputdirFig,filesep,'RF1.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF1.tif'])

h1 = figure('pos',POS);
ax1 = axes('parent',h1,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax1,'on');
% for i = 1:5
%     CirclePlot(i/5*maxpercent,ax1);
% end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)*maxpercent],[0,YShow(i,1)*maxpercent],'k--','parent',ax1);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
    for j = 1:5
        plot([XShow2(i,1)*j/5*maxpercent,XShow2(i+1,1)*j/5*maxpercent],...
            [YShow2(i,1)*j/5*maxpercent,YShow2(i+1,1)*j/5*maxpercent],...
            '--','parent',ax1,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)],...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)],...
        'parent',ax1,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i),YShow2(i,1)*SXpercent2(i),'parent',ax1,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax1,[-maxpercent,maxpercent,-maxpercent,maxpercent])
title('Perceont for Each Seed ROI')
saveas(h1,[OutputdirFig,filesep,'RF2.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF2.tif'])
%%

h2 = figure('pos',POS);
ax2 = axes('parent',h2,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax2,'on');
% for i = 1:5
%     CirclePlot(i/5*length(IX),ax2);
% end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)*length(IX)],[0,YShow(i,1)*length(IX)],'k--','parent',ax2);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
    for j = 1:5
        plot([XShow2(i,1)*j/5,XShow2(i+1,1)*j/5]*length(IX),...
            [YShow2(i,1)*j/5,YShow2(i+1,1)*j/5]*length(IX),...
            '--','parent',ax2,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        'parent',ax2,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i)*length(IX),YShow2(i,1)*SXpercent2(i)*length(IX),'parent',ax2,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax2,[-1*length(IX),1*length(IX),-1*length(IX),1*length(IX)])
title('Voxel Numbers for Each Seed ROI')

saveas(h2,[OutputdirFig,filesep,'RF3.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF3.tif'])

h3 = figure('pos',POS);
ax3 = axes('parent',h3,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax3,'on');
% for i = 1:5
%     CirclePlot(i/5*maxpercent*length(IX),ax3);
% end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)*maxpercent]*length(IX),[0,YShow(i,1)*maxpercent]*length(IX),'k--','parent',ax3);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
    for j = 1:5
        plot([XShow2(i,1)*j/5*maxpercent*length(IX),XShow2(i+1,1)*j/5*maxpercent*length(IX)],...
            [YShow2(i,1)*j/5*maxpercent*length(IX),YShow2(i+1,1)*j/5*maxpercent*length(IX)],...
            '--','parent',ax3,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        'parent',ax3,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i)*length(IX),YShow2(i,1)*SXpercent2(i)*length(IX),'parent',ax3,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax3,[-maxpercent,maxpercent,-maxpercent,maxpercent]*length(IX))
title('Voxel Numbers for Each Seed ROI')
saveas(h3,[OutputdirFig,filesep,'RF4.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF4.tif'])
%%

h = figure('pos',POS);
ax = axes('parent',h,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax,'on');
for i = 1:5
    CirclePlot(i/5,ax);
end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)],[0,YShow(i,1)],'k--','parent',ax);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
%     for j = 1:5
%         plot([XShow2(i,1)*j/5,XShow2(i+1,1)*j/5],...
%             [YShow2(i,1)*j/5,YShow2(i+1,1)*j/5],...
%             '--','parent',ax,'linewidth',0.5,'color',[0.5,0.5,0.5]);
%     end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)],...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)],...
        'parent',ax,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i),YShow2(i,1)*SXpercent2(i),'parent',ax,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax,[-1,1,-1,1])
title('Perceont for Each Seed ROI')

saveas(h,[OutputdirFig,filesep,'RF5.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF5.tif'])

h1 = figure('pos',POS);
ax1 = axes('parent',h1,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax1,'on');
for i = 1:5
    CirclePlot(i/5*maxpercent,ax1);
end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)*maxpercent],[0,YShow(i,1)*maxpercent],'k--','parent',ax1);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
%     for j = 1:5
%         plot([XShow2(i,1)*j/5*maxpercent,XShow2(i+1,1)*j/5*maxpercent],...
%             [YShow2(i,1)*j/5*maxpercent,YShow2(i+1,1)*j/5*maxpercent],...
%             '--','parent',ax1,'linewidth',0.5,'color',[0.5,0.5,0.5]);
%     end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)],...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)],...
        'parent',ax1,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i),YShow2(i,1)*SXpercent2(i),'parent',ax1,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax1,[-maxpercent,maxpercent,-maxpercent,maxpercent])
title('Perceont for Each Seed ROI')
saveas(h1,[OutputdirFig,filesep,'RF6.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF6.tif'])
%%

h2 = figure('pos',POS);
ax2 = axes('parent',h2,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax2,'on');
for i = 1:5
    CirclePlot(i/5*length(IX),ax2);
end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)*length(IX)],[0,YShow(i,1)*length(IX)],'k--','parent',ax2);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
%     for j = 1:5
%         plot([XShow2(i,1)*j/5,XShow2(i+1,1)*j/5]*length(IX),...
%             [YShow2(i,1)*j/5,YShow2(i+1,1)*j/5]*length(IX),...
%             '--','parent',ax2,'linewidth',0.5,'color',[0.5,0.5,0.5]);
%     end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        'parent',ax2,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i)*length(IX),YShow2(i,1)*SXpercent2(i)*length(IX),'parent',ax2,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax2,[-1*length(IX),1*length(IX),-1*length(IX),1*length(IX)])
title('Voxel Numbers for Each Seed ROI')

saveas(h2,[OutputdirFig,filesep,'RF7.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF7.tif'])

h3 = figure('pos',POS);
ax3 = axes('parent',h3,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax3,'on');
for i = 1:5
    CirclePlot(i/5*maxpercent*length(IX),ax3);
end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)*maxpercent]*length(IX),[0,YShow(i,1)*maxpercent]*length(IX),'k--','parent',ax3);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
%     for j = 1:5
%         plot([XShow2(i,1)*j/5*maxpercent*length(IX),XShow2(i+1,1)*j/5*maxpercent*length(IX)],...
%             [YShow2(i,1)*j/5*maxpercent*length(IX),YShow2(i+1,1)*j/5*maxpercent*length(IX)],...
%             '--','parent',ax3,'linewidth',0.5,'color',[0.5,0.5,0.5]);
%     end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        'parent',ax3,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i)*length(IX),YShow2(i,1)*SXpercent2(i)*length(IX),'parent',ax3,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax3,[-maxpercent,maxpercent,-maxpercent,maxpercent]*length(IX))
title('Voxel Numbers for Each Seed ROI')

saveas(h3,[OutputdirFig,filesep,'RF8.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF8.tif'])
%%

h = figure('pos',POS);
ax = axes('parent',h,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax,'on');
for i = 5
    CirclePlot(i/5,ax);
end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)],[0,YShow(i,1)],'k--','parent',ax);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
    for j = 1:5
        plot([XShow2(i,1)*j/5,XShow2(i+1,1)*j/5],...
            [YShow2(i,1)*j/5,YShow2(i+1,1)*j/5],...
            '--','parent',ax,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)],...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)],...
        'parent',ax,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i),YShow2(i,1)*SXpercent2(i),'parent',ax,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax,[-1,1,-1,1])
title('Perceont for Each Seed ROI')

saveas(h,[OutputdirFig,filesep,'RF9.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF9.tif'])

h1 = figure('pos',POS);
ax1 = axes('parent',h1,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax1,'on');
for i = 5
    CirclePlot(i/5*maxpercent,ax1);
end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)*maxpercent],[0,YShow(i,1)*maxpercent],'k--','parent',ax1);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
    for j = 1:5
        plot([XShow2(i,1)*j/5*maxpercent,XShow2(i+1,1)*j/5*maxpercent],...
            [YShow2(i,1)*j/5*maxpercent,YShow2(i+1,1)*j/5*maxpercent],...
            '--','parent',ax1,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)],...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)],...
        'parent',ax1,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i),YShow2(i,1)*SXpercent2(i),'parent',ax1,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax1,[-maxpercent,maxpercent,-maxpercent,maxpercent])
title('Perceont for Each Seed ROI')
saveas(h1,[OutputdirFig,filesep,'RF10.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF10.tif'])
%%

h2 = figure('pos',POS);
ax2 = axes('parent',h2,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax2,'on');
for i = 5
    CirclePlot(i/5*length(IX),ax2);
end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)*length(IX)],[0,YShow(i,1)*length(IX)],'k--','parent',ax2);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
    for j = 1:5
        plot([XShow2(i,1)*j/5,XShow2(i+1,1)*j/5]*length(IX),...
            [YShow2(i,1)*j/5,YShow2(i+1,1)*j/5]*length(IX),...
            '--','parent',ax2,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        'parent',ax2,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i)*length(IX),YShow2(i,1)*SXpercent2(i)*length(IX),'parent',ax2,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax2,[-1*length(IX),1*length(IX),-1*length(IX),1*length(IX)])
title('Voxel Numbers for Each Seed ROI')

saveas(h2,[OutputdirFig,filesep,'RF11.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF11.tif'])

h3 = figure('pos',POS);
ax3 = axes('parent',h3,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax3,'on');
for i = 5
    CirclePlot(i/5*maxpercent*length(IX),ax3);
end
for i = 1:seednum
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)*maxpercent]*length(IX),[0,YShow(i,1)*maxpercent]*length(IX),'k--','parent',ax3);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:seednum    
    for j = 1:5
        plot([XShow2(i,1)*j/5*maxpercent*length(IX),XShow2(i+1,1)*j/5*maxpercent*length(IX)],...
            [YShow2(i,1)*j/5*maxpercent*length(IX),YShow2(i+1,1)*j/5*maxpercent*length(IX)],...
            '--','parent',ax3,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)]*length(IX),...
        'parent',ax3,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i)*length(IX),YShow2(i,1)*SXpercent2(i)*length(IX),'parent',ax3,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax3,[-maxpercent,maxpercent,-maxpercent,maxpercent]*length(IX))
title('Voxel Numbers for Each Seed ROI')

saveas(h3,[OutputdirFig,filesep,'RF12.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'RF12.tif'])
%%
h = figure('pos',POS);

ax = axes('parent',h,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax,'on');
for i = 1:5
    AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)],[0,YShow(i,1)],'k--','parent',ax);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
SXpercent2 = [SXpercent;SXpercent(1)];
for i = 1:5    
    for j = 1:5
        plot([XShow2(i,1)*j/5,XShow2(i+1,1)*j/5],...
            [YShow2(i,1)*j/5,YShow2(i+1,1)*j/5],...
            '--','parent',ax,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
    plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)],...
        [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)],...
        'parent',ax,...
        'linewidth',2,...
        'color','k');
    plot(XShow2(i,1)*SXpercent2(i),YShow2(i,1)*SXpercent2(i),'parent',ax,...
        'marker','o',...
        'markersize',6,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:))
end
axis(ax,[-1,1,-1,1])
title('Percent for Each Seed ROI')

saveas(h3,[OutputdirFig,filesep,'Labs.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'Labs.tif'])

end
function CirclePlot(R,ax)
alpha=0:pi/20:2*pi;    %½Ç¶È[0,2*pi]

% R=2;                   %°ë¾¶
x=R*cos(alpha);
y=R*sin(alpha);
plot(x,y,'--','parent',ax,'linewidth',0.5,'color',[0.5,0.5,0.5]);
end