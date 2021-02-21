function PlotRadarGraph_multi(Data,seednum,maxpercent,Patchcolor,markersize,sepnum)
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
for isum = 1:size(Data,2);
    Datanew = Data(:,isum);
    patchcolor = Patchcolor(isum,:);
    SXpercent = Datanew/sum(Datanew);
    for i = 1:seednum
        AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
        XShow(i,1) = cos(AlphaVal(i,1));
        YShow(i,1) = sin(AlphaVal(i,1));
    end
    XShow2 = [XShow;XShow(1)];
    YShow2 = [YShow;YShow(1)];
    SXpercent2 = [SXpercent;SXpercent(1)];
    patch(XShow2.*SXpercent2,YShow2.*SXpercent2,patchcolor)
    for i = 1:seednum
        plot([0,XShow(i,1)*maxpercent],[0,YShow(i,1)*maxpercent],'k--','parent',ax);
    end
    for i = 1:seednum
        for j = 1:sepnum
            plot([XShow2(i,1)*j/sepnum*maxpercent,XShow2(i+1,1)*j/sepnum*maxpercent],...
                [YShow2(i,1)*j/sepnum*maxpercent,YShow2(i+1,1)*j/sepnum*maxpercent],...
                '--','parent',ax,'linewidth',1,'color',[0,0,0]);
        end
        plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)],...
            [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)],...
            'parent',ax,...
            'linewidth',2,...
            'color','k');
    end
    for i = 1:seednum
        plot(XShow2(i,1)*SXpercent2(i),YShow2(i,1)*SXpercent2(i),'parent',ax,...
            'marker','o',...
            'markersize',markersize,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
            'MarkerEdgeColor',SEEDCOLORSHOW(i,:),...
            'linewidth',2,'color','k',...
            'MarkerEdgeColor','k')
    end
    axis(ax,[-1,1,-1,1]*maxpercent*1.1)
end

h = figure('pos',POS);
ax = axes('parent',h,...
    'pos',[0.1,0.1,0.8,0.8]);
hold(ax,'on');
for isum = 1:size(Data,2);
    Datanew = Data(:,isum);
    patchcolor = Patchcolor(isum,:);
    SXpercent = Datanew/sum(Datanew);
    
    for i = 1:seednum
        AlphaVal(i,1) = pi/2-2*pi/seednum*(i-1);
        XShow(i,1) = cos(AlphaVal(i,1));
        YShow(i,1) = sin(AlphaVal(i,1));
    end
    XShow2 = [XShow;XShow(1)];
    YShow2 = [YShow;YShow(1)];
    SXpercent2 = [SXpercent;SXpercent(1)];
    patch(XShow2.*SXpercent2,YShow2.*SXpercent2,patchcolor)
    for i = 1:seednum
        plot([0,XShow(i,1)*maxpercent],[0,YShow(i,1)*maxpercent],'k--','parent',ax);
    end
    for i = 1:sepnum
        CirclePlot(i/sepnum*maxpercent,ax);
    end
    for i = 1:seednum
        plot([XShow2(i,1)*SXpercent2(i),XShow2(i+1,1)*SXpercent2(i+1)],...
            [YShow2(i,1)*SXpercent2(i),YShow2(i+1,1)*SXpercent2(i+1)],...
            'parent',ax,...
            'linewidth',2,...
            'color','k');
    end
    for i = 1:seednum
        plot(XShow2(i,1)*SXpercent2(i),YShow2(i,1)*SXpercent2(i),'parent',ax,...
            'marker','o',...
            'markersize',markersize,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
            'MarkerEdgeColor',SEEDCOLORSHOW(i,:),...
            'linewidth',2,'color','k',...
            'MarkerEdgeColor','k')
    end
    axis(ax,[-1,1,-1,1]*maxpercent*1.1)
end
end

function CirclePlot(R,ax)
alpha=0:pi/20:2*pi;    %½Ç¶È[0,2*pi]

% R=2;                   %°ë¾¶
x=R*cos(alpha);
y=R*sin(alpha);
plot(x,y,'--','parent',ax,'linewidth',1,'color',[0 0 0]);
end