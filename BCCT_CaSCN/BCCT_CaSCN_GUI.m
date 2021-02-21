function BCCT_CaSCN_GUI
CaSCN.fig = figure('unit','norm',...
    'pos',[0.4,0.4,0.3,0.2],...
    'menubar','none',...
    'Name','Causal network of Structural Covariance Connectivity(CaSCN)');

CaSCN.Mappb = uicontrol('unit','norm',...
    'pos',[0.125,0.3,0.15,0.4],'style','pushbutton','string','Map(Volume)');
CaSCN.Matpb = uicontrol('unit','norm',...
    'pos',[0.325,0.3,0.15,0.4],'style','pushbutton','string','Mat(Volume)');
CaSCN.Surfacepb = uicontrol('unit','norm',...
    'pos',[0.525,0.3,0.15,0.4],'style','pushbutton','string','Map(Surface)');
CaSCN.SurfaceROIpb = uicontrol('unit','norm',...
    'pos',[0.725,0.3,0.15,0.4],'style','pushbutton','string','Mat(Surface)');

CaSCN.return = uicontrol('unit','norm',...
    'pos',[0.7,0.1,0.1,0.15],'style','pushbutton','string','Return');
CaSCN.Exit = uicontrol('unit','norm',...
    'pos',[0.85,0.1,0.1,0.15],'style','pushbutton','string','Exit');

set(CaSCN.Mappb,'callback',{@SCNMAP,CaSCN});
set(CaSCN.Matpb,'callback',{@SCNMAT,CaSCN});
set(CaSCN.Surfacepb,'callback',{@SURFMAP,CaSCN});
set(CaSCN.SurfaceROIpb,'callback',{@SURFROI,CaSCN});
set(CaSCN.return,'callback',{@Return,CaSCN});
set(CaSCN.Exit,'callback',{@Exit,CaSCN});
end
function SCNMAP(varargin)
CaSCN = varargin{3};
close(CaSCN.fig);
BCCT_CaSCN_MAP_GUI;
end
function SCNMAT(varargin)
CaSCN = varargin{3};
close(CaSCN.fig);
BCCT_CaSCN_Mat_GUI;
end
function SURFMAP(varargin)
CaSCN = varargin{3};
close(CaSCN.fig);
BCCT_CaSCN_Surf_GUI;
end
function SURFROI(varargin)
CaSCN = varargin{3};
close(CaSCN.fig);
BCCT_CaSCN_SurfROI_GUI;
end
function Return(varargin)
CaSCN = varargin{3};
close(CaSCN.fig);
BCCT;
end
function Exit(varargin)
CaSCN = varargin{3};
close(CaSCN.fig);
end