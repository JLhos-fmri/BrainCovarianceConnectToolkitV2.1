function BCCT_Modulate_GUI

MOD.fig = figure('unit','norm',...
    'pos',[0.4,0.4,0.3,0.2],...
    'menubar','none',...
    'Name','Modulation on Covariance Connectivity');

MOD.Mappb = uicontrol('unit','norm',...
    'pos',[0.125,0.3,0.15,0.4],'style','pushbutton','string','Map(Volume)');
MOD.Matpb = uicontrol('unit','norm',...
    'pos',[0.325,0.3,0.15,0.4],'style','pushbutton','string','Mat(Volume)');
MOD.Surfacepb = uicontrol('unit','norm',...
    'pos',[0.525,0.3,0.15,0.4],'style','pushbutton','string','Map(Surface)');
MOD.SurfaceROIpb = uicontrol('unit','norm',...
    'pos',[0.725,0.3,0.15,0.4],'style','pushbutton','string','Mat(Surface)');

MOD.return = uicontrol('unit','norm',...
    'pos',[0.7,0.1,0.1,0.15],'style','pushbutton','string','Return');
MOD.Exit = uicontrol('unit','norm',...
    'pos',[0.85,0.1,0.1,0.15],'style','pushbutton','string','Exit');

set(MOD.Mappb,'callback',{@SCNMAP,MOD});
set(MOD.Matpb,'callback',{@SCNMAT,MOD});
set(MOD.Surfacepb,'callback',{@SURFMAP,MOD});
set(MOD.SurfaceROIpb,'callback',{@SURFROI,MOD});
set(MOD.return,'callback',{@Return,MOD});
set(MOD.Exit,'callback',{@Exit,MOD});
end
function SCNMAP(varargin)
MOD = varargin{3};
close(MOD.fig);
BCCT_MOD_MAP_GUI;
end
function SCNMAT(varargin)
MOD = varargin{3};
close(MOD.fig);
BCCT_MOD_Mat_GUI;
end
function SURFMAP(varargin)
MOD = varargin{3};
close(MOD.fig);
BCCT_MOD_Surf_GUI;
end
function SURFROI(varargin)
MOD = varargin{3};
close(MOD.fig);
BCCT_MOD_SurfROI_GUI;
end
function Return(varargin)
MOD = varargin{3};
close(MOD.fig);
BCCT;
end
function Exit(varargin)
MOD = varargin{3};
close(MOD.fig);
end

