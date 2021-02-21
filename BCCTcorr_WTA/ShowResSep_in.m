function ShowResSep_in(Outputdir)
% Outputdir = uigetdir(pwd,'Output Directory Selection');
pathfile = which('ASBCmain.m');
[pat nam ext] = fileparts(pathfile);
BackGroundFile = fullfile(pat,'SomeTemplates','mni_icbm152_t1.nii');
[vBackGroud datBackGroud] = Dynamic_read_dir_NIFTI(BackGroundFile);
OutputdirFig = [Outputdir,filesep,'ResFigure'];
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
pathmfile = which('AS_WTA_GUI.m');
[pat nam ext] = fileparts(pathmfile);
load(fullfile(pat,'SEEDCOLOR.mat'))
load(fullfile(pat,'JETcode.mat'));
%%
h0 = figure;
imagesc(r',[-0.5,0.5]);
colormap(JETcode);
colorbar
title('Orig Rmap');
saveas(h0,[OutputdirFig,filesep,'OrigRmap1.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'OrigRmap1.tif'])

[IX IY] = sort(maxind);
LABnum = seednum;
sep = 2/(LABnum+1);
valshow = 1:sep:3;
H1 = figure;
R_SHOW = r(:,IY)';
R_SHOW2 = R_SHOW;
R_SHOW2(R_SHOW<-0.5) = -0.5;
R_SHOW2(R_SHOW>0.5) = 0.5;
NUMSEP1 = size(JETcode,1)*2;
K = floor(NUMSEP1/seednum);
for i = 1:seednum
    FORSHOWCOLORBAR(1+K*(i-1):K*i,1) = SEEDCOLORSHOW(i,1);
    FORSHOWCOLORBAR(1+K*(i-1):K*i,2) = SEEDCOLORSHOW(i,2);
    FORSHOWCOLORBAR(1+K*(i-1):K*i,3) = SEEDCOLORSHOW(i,3);
end
FORSHOWCOLORLow = zeros(128,3);
FORSHOWCOLORLow(33:32*3,:) = JETcode;
FORSHOWCOLORLow(1:31,1) = JETcode(1,1);
FORSHOWCOLORLow(1:31,2) = JETcode(1,2);
FORSHOWCOLORLow(1:31,3) = JETcode(1,3);
FORSHOWCOLORLow(97:128,1) = JETcode(64,1);
FORSHOWCOLORLow(97:128,2) = JETcode(64,2);
FORSHOWCOLORLow(97:128,3) = JETcode(64,3);
FOWSHOWCOLORfinal = [FORSHOWCOLORLow;FORSHOWCOLORBAR];
IXs = IX';
IXforshow = zeros(size(IXs));
for i = 1:seednum
    IXforshow(IXs==i) = valshow(i+1);
end
Rshow = [IXforshow,R_SHOW2];
imagesc(Rshow,[-1,3]);
colormap(FOWSHOWCOLORfinal);
colorbar;
hold on;
for i = 1:seednum-1
    lineind = find(IX==i+1);
    if ~isempty(lineind)
        linepos = lineind(1)-0.5;
        plot([0.5,seednum+1.5],[linepos,linepos],'linewidth',2,'color',[0 0 0]);
    end
end
title('Sorted Rmap1');
saveas(H1,[OutputdirFig,filesep,'SortedRmap1.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'SortedRmap1.tif'])
%%
H2 = figure;
R_SHOW = r(:,IY)';
R_SHOW2 = R_SHOW;
NUMSEP1 = size(JETcode,1);
K = floor(NUMSEP1/seednum);
clear FORSHOWCOLORBAR
for i = 1:seednum
    FORSHOWCOLORBAR(1+K*(i-1):K*i,1) = SEEDCOLORSHOW(i,1);
    FORSHOWCOLORBAR(1+K*(i-1):K*i,2) = SEEDCOLORSHOW(i,2);
    FORSHOWCOLORBAR(1+K*(i-1):K*i,3) = SEEDCOLORSHOW(i,3);
end
FORSHOWCOLORLow = JETcode;
FOWSHOWCOLORfinal = [FORSHOWCOLORLow;FORSHOWCOLORBAR];
IXs = IX';
IXforshow = zeros(size(IXs));
for i = 1:seednum
    IXforshow(IXs==i) = valshow(i+1);
end
Rshow = [IXforshow,R_SHOW2];
imagesc(Rshow,[-1,3]);
colormap(FOWSHOWCOLORfinal);
colorbar;
hold on;
for i = 1:seednum-1
    lineind = find(IX==i+1);
    if ~isempty(lineind)
        linepos = lineind(1)-0.5;
        plot([0.5,seednum+1.5],[linepos,linepos],'linewidth',2,'color',[0 0 0]);
    end
end
title('Sorted Rmap2');
saveas(H2,[OutputdirFig,filesep,'SortedRmap2.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'SortedRmap2.tif'])
%% slice view
[vtarget vdat] = Dynamic_read_dir_NIFTI(RealCompPara.Targetdir);
outnameid = fullfile(Outputdir,'mixedMaxID.nii');
[vout voutdat] = Dynamic_read_dir_NIFTI(outnameid);
indout = find(voutdat);
VOUTSHOW = reshape(voutdat,vout.dim);
voutindexexist = unique(voutdat);
voutindexexist(1) = [];

[indx indy indz] = ind2sub(vout.dim,indout);
indxrange = [min(indx),max(indx)];
indyrange = [min(indy),max(indy)];
indzrange = [min(indz),max(indz)];

numshow = 15;
indxind = indxrange(1):(indxrange(2)-indxrange(1))/(numshow-1):indxrange(2);
indxind = floor(indxind);
indyind = indyrange(1):(indyrange(2)-indyrange(1))/(numshow-1):indyrange(2);
indyind = floor(indyind);
indzind = indzrange(1):(indzrange(2)-indzrange(1))/(numshow-1):indzrange(2);
indzind = floor(indzind);
startmni = cor2mni([indxrange(1),indyrange(1),indzrange(1)],vout.mat);
endmin = cor2mni([indxrange(2),indyrange(2),indzrange(2)],vout.mat);
xind = floor(startmni(1)):(floor(endmin(1))-floor(startmni(1)))/(numshow-1):floor(endmin(1));
xind = floor(xind);
yind = floor(startmni(2)):(floor(endmin(2))-floor(startmni(2)))/(numshow-1):floor(endmin(2));
yind = floor(yind);
zind = floor(startmni(3)):(floor(endmin(3))-floor(startmni(3)))/(numshow-1):floor(endmin(3));
zind = floor(zind);

Sview.fig = figure('Name','Slice Viewer',...       
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'color',[0.95 0.95 0.95],...
    'position',[0.05 0.05 0.9 0.9]);
Sview.uibutton(1) = uibuttongroup('parent',Sview.fig,...
    'units','normalized',...
    'pos',[0.05,0.05,0.9,0.8]);
Sview.uibutton(2) = uibuttongroup('parent',Sview.fig,...
    'units','normalized',...
    'pos',[0.05,0.05,0.9,0.8]);
Sview.uibutton(3) = uibuttongroup('parent',Sview.fig,...
    'units','normalized',...
    'pos',[0.05,0.05,0.9,0.8]);
Sview.viewopt(1) = uicontrol('parent',Sview.fig,...
    'units','normalized',...
    'pos',[0.05,0.9,0.2,0.08],...
    'style','radiobutton',...
    'string','Axial',...
    'fontunits', 'normalized',...
    'fontsize',0.6,...
    'fontweight','bold',...
    'horizontalalign','center',...
    'value',1);
Sview.viewopt(2) = uicontrol('parent',Sview.fig,...
    'units','normalized',...
    'pos',[0.3,0.9,0.2,0.08],...
    'style','radiobutton',...
    'string','Coronal',...
    'fontunits', 'normalized',...
    'fontsize',0.6,...
    'fontweight','bold',...
    'horizontalalign','center',...
    'value',0);
Sview.viewopt(3) = uicontrol('parent',Sview.fig,...
    'units','normalized',...
    'pos',[0.6,0.9,0.2,0.08],...
    'style','radiobutton',...
    'string','Sagittal',...
    'fontunits', 'normalized',...
    'fontsize',0.6,...
    'fontweight','bold',...
    'horizontalalign','center',...
    'value',0);
SepPart = seednum+1;
colormapshow = [0 0 0;SEEDCOLORSHOW(1:seednum,:);1 1 1];
for i = 1:numshow
    Sview.lab(i) = uicontrol('parent',Sview.uibutton(1),...
        'units','norm',...
        'pos',[0.05+(0.9/numshow)*(i-1),0.85,0.9/numshow,0.05],...
        'style','text',...
        'string',['z = ',num2str(zind(i))],...
        'fontunits', 'normalized',...
        'fontsize',0.6,...
        'fontweight','bold',...
        'horizontalalign','center');
    for j = 1:SepPart
        Sview.SAxes(i,j) = axes('parent',Sview.uibutton(1),...
            'pos',[0.05+(0.9/numshow)*(i-1),0.85-0.8/SepPart*j,0.9/numshow,0.8/SepPart]);
%             'pos',[0.05+0.09*(i-1),0.05+0.8/SepPart*(j-1),0.09,0.8/SepPart]);
        axis(Sview.SAxes(i,j),'off');
    end
    DATshow = rot90(squeeze(VOUTSHOW(indxrange(1)-1:indxrange(2)+1,indyrange(1)-1:indyrange(2)+1,indzind(i))));
    imshow(DATshow,[0,6],'parent',Sview.SAxes(i,1));
    axis(Sview.SAxes(i,1),[0.5,size(DATshow,2)+0.5,0.5,size(DATshow,1)+0.5])
    hold(Sview.SAxes(i,1),'on');
    plot([0.5,0.5],[0.5,size(DATshow,1)+0.5],'parent',Sview.SAxes(i,1),'linewidth',2,'color',[1 1 1])
    plot([size(DATshow,2)+0.5,size(DATshow,2)+0.5],[0.5,size(DATshow,1)+0.5],'parent',Sview.SAxes(i,1),'linewidth',2,'color',[1 1 1])
    plot([0.5,size(DATshow,2)+0.5],[0.5,0.5],'parent',Sview.SAxes(i,1),'linewidth',2,'color',[1 1 1])
    plot([0.5,size(DATshow,2)+0.5],[size(DATshow,1)+0.5,size(DATshow,1)+0.5],'parent',Sview.SAxes(i,1),'linewidth',2,'color',[1 1 1])
    hold(Sview.SAxes(i,1),'off');
    for j = 1:seednum
        DATshow2 = ones(size(DATshow))*(seednum+1);
        DATshow2(DATshow==voutindexexist(j)) = j;
        DATshow2(DATshow==0) = 0;
        imshow(DATshow2,[0,6],'parent',Sview.SAxes(i,j+1));
        axis(Sview.SAxes(i,j+1),[0.5,size(DATshow,2)+0.5,0.5,size(DATshow,1)+0.5])
        
        hold(Sview.SAxes(i,j+1),'on');
        plot([0.5,0.5],[0.5,size(DATshow,1)+0.5],'parent',Sview.SAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        plot([size(DATshow,2)+0.5,size(DATshow,2)+0.5],[0.5,size(DATshow,1)+0.5],'parent',Sview.SAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        plot([0.5,size(DATshow,2)+0.5],[0.5,0.5],'parent',Sview.SAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        plot([0.5,size(DATshow,2)+0.5],[size(DATshow,1)+0.5,size(DATshow,1)+0.5],'parent',Sview.SAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        hold(Sview.SAxes(i,j+1),'off');
    end
end
colormap(colormapshow)
%%
for i = 1:numshow
    Sview.Ylab(i) = uicontrol('parent',Sview.uibutton(2),...
        'units','norm',...
        'pos',[0.05+(0.9/numshow)*(i-1),0.85,0.9/numshow,0.05],...
        'style','text',...
        'string',['y = ',num2str(yind(i))],...
        'fontunits', 'normalized',...
        'fontsize',0.6,...
        'fontweight','bold',...
        'horizontalalign','center');
    for j = 1:SepPart
        Sview.YSAxes(i,j) = axes('parent',Sview.uibutton(2),...
            'pos',[0.05+(0.9/numshow)*(i-1),0.85-0.8/SepPart*j,0.9/numshow,0.8/SepPart]);
%             'pos',[0.05+0.09*(i-1),0.05+0.8/SepPart*(j-1),0.09,0.8/SepPart]);
        axis(Sview.YSAxes(i,j),'off');
    end
    DATshow = rot90(squeeze(VOUTSHOW(indxrange(1)-1:indxrange(2)+1,indyind(i),indzrange(1)-1:indzrange(2)+1)));
    imshow(DATshow,[0,6],'parent',Sview.YSAxes(i,1));
%     AA = get(Sview.YSAxes(i,1),'xlim')
    axis(Sview.YSAxes(i,1),[0.5,size(DATshow,2)+0.5,0.5,size(DATshow,1)+0.5])
    hold(Sview.YSAxes(i,1),'on');
    plot([0.5,0.5],[0.5,size(DATshow,1)+0.5],'parent',Sview.YSAxes(i,1),'linewidth',2,'color',[1 1 1])
    plot([size(DATshow,2)+0.5,size(DATshow,2)+0.5],[0.5,size(DATshow,1)+0.5],'parent',Sview.YSAxes(i,1),'linewidth',2,'color',[1 1 1])
    plot([0.5,size(DATshow,2)+0.5],[0.5,0.5],'parent',Sview.YSAxes(i,1),'linewidth',2,'color',[1 1 1])
    plot([0.5,size(DATshow,2)+0.5],[size(DATshow,2)+0.5,size(DATshow,1)+0.5],'parent',Sview.YSAxes(i,1),'linewidth',2,'color',[1 1 1])
    for j = 1:seednum
        DATshow2 = ones(size(DATshow))*(seednum+1);
        DATshow2(DATshow==voutindexexist(j)) = j;
        DATshow2(DATshow==0) = 0;
        imshow(DATshow2,[0,6],'parent',Sview.YSAxes(i,j+1));
        axis(Sview.YSAxes(i,j+1),[0.5,size(DATshow,2)+0.5,0.5,size(DATshow,1)+0.5])
        
        hold(Sview.YSAxes(i,j+1),'on');
        plot([0.5,0.5],[0.5,size(DATshow,1)+0.5],'parent',Sview.YSAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        plot([size(DATshow,2)+0.5,size(DATshow,2)+0.5],[0.5,size(DATshow,1)+0.5],'parent',Sview.YSAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        plot([0.5,size(DATshow,2)+0.5],[0.5,0.5],'parent',Sview.YSAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        plot([0.5,size(DATshow,2)+0.5],[size(DATshow,1)+0.5,size(DATshow,1)+0.5],'parent',Sview.YSAxes(i,j+1),'linewidth',2,'color',[1 1 1])
    end
end
colormap(colormapshow)
set(Sview.uibutton(2),'vis','off')
%%
for i = 1:numshow
    Sview.Xlab(i) = uicontrol('parent',Sview.uibutton(3),...
        'units','norm',...
        'pos',[0.05+(0.9/numshow)*(i-1),0.85,0.9/numshow,0.05],...
        'style','text',...
        'string',['x = ',num2str(xind(i))],...
        'fontunits', 'normalized',...
        'fontsize',0.6,...
        'fontweight','bold',...
        'horizontalalign','center');
    for j = 1:SepPart
        Sview.XSAxes(i,j) = axes('parent',Sview.uibutton(3),...
            'pos',[0.05+(0.9/numshow)*(i-1),0.85-0.8/SepPart*j,0.9/numshow,0.8/SepPart]);
%             'pos',[0.05+0.09*(i-1),0.05+0.8/SepPart*(j-1),0.09,0.8/SepPart]);
        axis(Sview.XSAxes(i,j),'off');
    end
    DATshowX = rot90(squeeze(VOUTSHOW(indxind(i),indyrange(1)-1:indyrange(2)+1,indzrange(1)-1:indzrange(2)+1)));
    imshow(DATshowX,[0,6],'parent',Sview.XSAxes(i,1));
    axis(Sview.XSAxes(i,1),[0.5,size(DATshowX,2)+0.5,0.5,size(DATshowX,1)+0.5])
    hold(Sview.XSAxes(i,1),'on');
    plot([0.5,0.5],[0.5,size(DATshowX,1)+0.5],'parent',Sview.XSAxes(i,1),'linewidth',2,'color',[1 1 1])
    plot([size(DATshowX,2)+0.5,size(DATshowX,2)+0.5],[0.5,size(DATshowX,1)+0.5],'parent',Sview.XSAxes(i,1),'linewidth',2,'color',[1 1 1])
    plot([0.5,size(DATshowX,2)+0.5],[0.5,0.5],'parent',Sview.XSAxes(i,1),'linewidth',2,'color',[1 1 1])
    plot([0.5,size(DATshowX,2)+0.5],[size(DATshowX,1)+0.5,size(DATshowX,1)+0.5],'parent',Sview.XSAxes(i,1),'linewidth',2,'color',[1 1 1])
    for j = 1:seednum
        DATshow2 = ones(size(DATshowX))*(seednum+1);
        DATshow2(DATshowX==voutindexexist(j)) = j;
        DATshow2(DATshowX==0) = 0;
        imshow(DATshow2,[0,6],'parent',Sview.XSAxes(i,j+1));
        axis(Sview.XSAxes(i,j+1),[0.5,size(DATshowX,2)+0.5,0.5,size(DATshowX,1)+0.5])
        
        hold(Sview.XSAxes(i,j+1),'on');
        plot([0.5,0.5],[0.5,size(DATshowX,1)+0.5],'parent',Sview.XSAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        plot([size(DATshowX,2)+0.5,size(DATshowX,2)+0.5],[0.5,size(DATshowX,1)+0.5],'parent',Sview.XSAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        plot([0.5,size(DATshowX,2)+0.5],[0.5,0.5],'parent',Sview.XSAxes(i,j+1),'linewidth',2,'color',[1 1 1])
        plot([0.5,size(DATshowX,2)+0.5],[size(DATshowX,1)+0.5,size(DATshowX,1)+0.5],'parent',Sview.XSAxes(i,j+1),'linewidth',2,'color',[1 1 1])
    end
end
colormap(colormapshow)
set(Sview.uibutton(3),'vis','off')
Sview.voutindexexist = voutindexexist;
%%
showfig = figure;
ax1 = axes('parent',showfig,...
    'pos',[0.05,0.05,0.25,0.9]);
ax2 = axes('parent',showfig,...
    'pos',[0.35,0.05,0.6,0.9]);
hold(ax1,'on');
hold(ax2,'on');
for i = 1:seednum
    SX(i,1) = length(find(IX==i));
end
SY = cumsum(SX);
COLORmark = SEEDCOLORSHOW(seednum:-1:1,:);
for i = 1:seednum
    bar(1,SY(seednum-i+1),'barwidth',0.2,...
        'FaceColor',COLORmark(i,:),...
        'parent',ax1);
    bar(i,SX(i),'barwidth',0.4,...
        'FaceColor',SEEDCOLORSHOW(i,:),...
        'parent',ax2)
end

saveas(showfig,[OutputdirFig,filesep,'voxelnums2.fig']);
print('-dtiff','-r300',[OutputdirFig,filesep,'voxelnums2.tif'])
%%
set(Sview.viewopt(1),'callback',{@AxialView,Sview,numshow,SepPart})
set(Sview.viewopt(2),'callback',{@CoronalView,Sview,numshow,SepPart})
set(Sview.viewopt(3),'callback',{@SagittalView,Sview,numshow,SepPart})
end

function AxialView(varargin)
Sview = varargin{3};
numshow = varargin{4};
SepPart = varargin{5};
voutindexexist = Sview.voutindexexist;
%%
set(Sview.viewopt(1),'val',1);
set(Sview.viewopt(2),'val',0);
set(Sview.viewopt(3),'val',0);
% for i = 1:numshow
%     set(Sview.Xlab(i),'vis','off');
%     for j = 1:SepPart
%         set(Sview.XSAxes(i,j),'vis','off')
%     end
% end
% for i = 1:numshow
%     set(Sview.Ylab(i),'vis','off');
%     for j = 1:SepPart
%         set(Sview.YSAxes(i,j),'vis','off')
%     end
% end
% 
% for i = 1:numshow
%     set(Sview.lab(i),'vis','on');
%     for j = 1:SepPart
%         set(Sview.SAxes(i,j),'vis','on')
%     end
% end

set(Sview.uibutton(1),'vis','on')
set(Sview.uibutton(2),'vis','off')
set(Sview.uibutton(3),'vis','off')
end
function CoronalView(varargin)
Sview = varargin{3};
numshow = varargin{4};
SepPart = varargin{5};
voutindexexist = Sview.voutindexexist;
%%
set(Sview.viewopt(1),'val',0);
set(Sview.viewopt(2),'val',1);
set(Sview.viewopt(3),'val',0);

% for i = 1:numshow
%     set(Sview.Xlab(i),'vis','off');
%     for j = 1:SepPart
%         set(Sview.XSAxes(i,j),'vis','off')
%     end
% end
% for i = 1:numshow
%     set(Sview.Ylab(i),'vis','on');
%     for j = 1:SepPart
%         set(Sview.YSAxes(i,j),'vis','on')
%     end
% end
% 
% for i = 1:numshow
%     set(Sview.lab(i),'vis','off');
%     for j = 1:SepPart
%         set(Sview.SAxes(i,j),'vis','off')
%     end
% end
set(Sview.uibutton(1),'vis','off')
set(Sview.uibutton(2),'vis','on')
set(Sview.uibutton(3),'vis','off')
end
function SagittalView(varargin)
Sview = varargin{3};
numshow = varargin{4};
SepPart = varargin{5};
voutindexexist = Sview.voutindexexist;
%%
set(Sview.viewopt(1),'val',0);
set(Sview.viewopt(2),'val',0);
set(Sview.viewopt(3),'val',1);
% for i = 1:numshow
%     set(Sview.Xlab(i),'vis','on');
%     for j = 1:SepPart
%         set(Sview.XSAxes(i,j),'vis','on')
%     end
% end
% for i = 1:numshow
%     set(Sview.Ylab(i),'vis','off');
%     for j = 1:SepPart
%         set(Sview.YSAxes(i,j),'vis','off')
%     end
% end
% 
% for i = 1:numshow
%     set(Sview.lab(i),'vis','off');
%     for j = 1:SepPart
%         set(Sview.SAxes(i,j),'vis','off')
%     end
% end
set(Sview.uibutton(1),'vis','off')
set(Sview.uibutton(2),'vis','off')
set(Sview.uibutton(3),'vis','on')
end
