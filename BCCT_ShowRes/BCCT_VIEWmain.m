function BCCT_VIEWmain
Hsize = get(0,'screensize');
msize = min(Hsize(3:4))*0.8;
Hasview.fig = figure('pos',[Hsize(3)/2-msize/4,Hsize(4)/2-msize/8,msize/2,msize/4],'name','BCCT Viewer V1.0');
Hasview.FC = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.025 0.6 0.3 0.3],'style','pushbutton','string','Map(SCN&Modulate)');
Hasview.GC = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.025 0.1 0.3 0.3],'style','pushbutton','string','Map(GCA)');
Hasview.MatrixFC = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.35 0.6 0.3 0.3],'style','pushbutton','string','Matrix(SCN&Modulate)');
Hasview.MatrixGC = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.35 0.1 0.3 0.3],'style','pushbutton','string','Matrix(GCA)');
Hasview.WTA = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.675 0.6 0.3 0.3],'style','pushbutton','string','Winner-Take-All');
Hasview.Return = uicontrol('parent',Hasview.fig,'units','norm','pos',[0.675 0.1 0.3 0.3],'style','pushbutton','string','Return');
set(Hasview.FC,'callback',{@FCmapshow,Hasview})
set(Hasview.GC,'callback',{@GCmapshow,Hasview})
set(Hasview.MatrixFC,'callback',{@FCmatrixshow,Hasview})
set(Hasview.MatrixGC,'callback',{@GCmatrixshow,Hasview})
set(Hasview.WTA,'callback',{@WTAshow,Hasview})
set(Hasview.Return,'callback',{@Return,Hasview})
end
function FCmapshow(varargin)
Hasview = varargin{3};
close(Hasview.fig);
BCCT_ShowMap_GUI;
end
function GCmapshow(varargin)
Hasview = varargin{3};
% uiwait(msgbox('comming soon'));
close(Hasview.fig);
answ1 = questdlg('which type of GC?','which type of GC?','Residual','Coefficient','Residual');
if strcmp(answ1,'Residual')
    BCCT_ShowGCMap_GUI_res;
else
    BCCT_ShowGCMap_GUI_coef;
end
end
function FCmatrixshow(varargin)
Hasview = varargin{3};
close(Hasview.fig);
BCCT_ShowFCmatrixGUI
% BCCT_ShowMatrix_GUI;
end
function GCmatrixshow(varargin)
Hasview = varargin{3};
% uiwait(msgbox('comming soon'));
close(Hasview.fig);
BCCT_ShowGCmatrix_GUI;
end
function WTAshow(varargin)
Hasview = varargin{3};
close(Hasview.fig);
BCCT_ShowWTAGUI;
end
function Return(varargin)
Hasview = varargin{3};
close(Hasview.fig);
BCCT;
end