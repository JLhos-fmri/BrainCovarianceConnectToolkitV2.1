function BCCT_SCN_GUI
SCN.fig = figure('unit','norm',...
    'pos',[0.4,0.4,0.3,0.2],...
    'menubar','none',...
    'Name','Structural Covariance Network(SCN)');

SCN.Mappb = uicontrol('unit','norm',...
    'pos',[0.125,0.3,0.15,0.4],'style','pushbutton','string','Map(Volume)');
SCN.Matpb = uicontrol('unit','norm',...
    'pos',[0.325,0.3,0.15,0.4],'style','pushbutton','string','Mat(Volume)');
SCN.Surfacepb = uicontrol('unit','norm',...
    'pos',[0.525,0.3,0.15,0.4],'style','pushbutton','string','Map(Surface)');
SCN.SurfaceROIpb = uicontrol('unit','norm',...
    'pos',[0.725,0.3,0.15,0.4],'style','pushbutton','string','Mat(Surface)');

SCN.return = uicontrol('unit','norm',...
    'pos',[0.7,0.1,0.1,0.15],'style','pushbutton','string','Return');
SCN.Exit = uicontrol('unit','norm',...
    'pos',[0.85,0.1,0.1,0.15],'style','pushbutton','string','Exit');

set(SCN.Mappb,'callback',{@SCNMAP,SCN});
set(SCN.Matpb,'callback',{@SCNMAT,SCN});
set(SCN.Surfacepb,'callback',{@SURFMAP,SCN});
set(SCN.SurfaceROIpb,'callback',{@SURFROI,SCN});
set(SCN.return,'callback',{@Return,SCN});
set(SCN.Exit,'callback',{@Exit,SCN});
end
function SCNMAP(varargin)
SCN = varargin{3};
close(SCN.fig);
BCCT_CON_MAP_GUI;
end
function SCNMAT(varargin)
SCN = varargin{3};
close(SCN.fig);
BCCT_CON_Mat_GUI;
end
function SURFMAP(varargin)
SCN = varargin{3};
close(SCN.fig);
BCCT_CON_Surf_GUI;
end
function SURFROI(varargin)
SCN = varargin{3};
close(SCN.fig);
BCCT_CON_SurfROI_GUI;
end
function Return(varargin)
SCN = varargin{3};
close(SCN.fig);
BCCT;
end
function Exit(varargin)
SCN = varargin{3};
close(SCN.fig);
end