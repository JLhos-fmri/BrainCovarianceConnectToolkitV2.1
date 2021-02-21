function AS_ShowGCmatrixGUI
Hsize = get(0,'screensize');
MIDPOINT = [Hsize(3)/2,Hsize(4)/2];
Asize = [100*3,100+40];
MaxSIZE = [Hsize(3) Hsize(4)]*0.4;
factor = MaxSIZE./Asize;
factornew = min(factor);
POSSIZE = Asize*factornew;
SelMod.fig = figure('position',[MIDPOINT(1)-POSSIZE(1)/4,MIDPOINT(2)-POSSIZE(2)/4,POSSIZE(1)/2,POSSIZE(2)/2],'name','GC Matrix Show Mod selection');
SelMod.GCtype = uibuttongroup('parent',SelMod.fig,'units','norm',...
    'Pos',[0.1,0.65,0.8,0.25]);
SelMod.Pval = uibuttongroup('parent',SelMod.fig,'units','norm',...
    'Pos',[0.1,0.35,0.8,0.25]);
SelMod.GCtype1 = uicontrol('parent',SelMod.GCtype,'units','norm',...
    'pos',[0.1,0.1,0.35,0.8],...
    'style','rad',...
    'string','residual',...
    'val',1);
SelMod.GCtype2 = uicontrol('parent',SelMod.GCtype,'units','norm',...
    'pos',[0.55,0.1,0.35,0.8],...
    'style','rad',...
    'string','coefficient',...
    'val',0);
SelMod.Ptype1 = uicontrol('parent',SelMod.Pval,'units','norm',...
    'pos',[0.1,0.1,0.25,0.8],...
    'style','rad',...
    'string','NoPval',...
    'val',1);
SelMod.Ptype2 = uicontrol('parent',SelMod.Pval,'units','norm',...
    'pos',[0.4,0.1,0.25,0.8],...
    'style','rad',...
    'string','Pval',...
    'val',0);
SelMod.Ptype3 = uicontrol('parent',SelMod.Pval,'units','norm',...
    'pos',[0.7,0.1,0.25,0.8],...
    'style','rad',...
    'string','Permutation',...
    'val',0);

SelMod.SetUp = uicontrol('parent',SelMod.fig,'units','norm',...
    'pos',[0.1,0.1,0.3,0.2],...
    'style','pushbutton',...
    'string','Finish Selection!');
SelMod.Return = uicontrol('parent',SelMod.fig,'units','norm',...
    'pos',[0.6,0.1,0.3,0.2],...
    'style','pushbutton',...
    'string','Return!');

set(SelMod.GCtype1,'callback',{@GCTYPE1,SelMod});
set(SelMod.GCtype2,'callback',{@GCTYPE2,SelMod});
set(SelMod.Ptype1,'callback',{@Ptype1,SelMod});
set(SelMod.Ptype2,'callback',{@Ptype2,SelMod});
set(SelMod.Ptype3,'callback',{@Ptype3,SelMod});
set(SelMod.SetUp,'callback',{@SetUpfin,SelMod});
set(SelMod.Return,'callback',{@Return,SelMod});
end

function GCTYPE1(varargin)
SelMod = varargin{3};
set(SelMod.GCtype1,'val',1);
set(SelMod.GCtype2,'val',0);
set(SelMod.Ptype2,'enable','on');
end
function GCTYPE2(varargin)
SelMod = varargin{3};
uiwait(msgbox('Current only support GCorder=1'));
set(SelMod.GCtype2,'val',1);
set(SelMod.GCtype1,'val',0);
set(SelMod.Ptype2,'enable','off');
set(SelMod.Ptype1,'val',1);
set(SelMod.Ptype2,'val',0);
set(SelMod.Ptype3,'val',0);
end
function Ptype1(varargin)
SelMod = varargin{3};
set(SelMod.Ptype1,'val',1);
set(SelMod.Ptype2,'val',0);
set(SelMod.Ptype3,'val',0);
end
function Ptype2(varargin)
SelMod = varargin{3};
set(SelMod.Ptype1,'val',0);
set(SelMod.Ptype2,'val',1);
set(SelMod.Ptype3,'val',0);
end
function Ptype3(varargin)
SelMod = varargin{3};
set(SelMod.Ptype1,'val',0);
set(SelMod.Ptype2,'val',0);
set(SelMod.Ptype3,'val',1);
end
function SetUpfin(varargin)
SelMod = varargin{3};
GCtype = get(SelMod.GCtype1,'val');
Ptype1 = get(SelMod.Ptype1,'val');
Ptype2 = get(SelMod.Ptype2,'val');
Ptype3 = get(SelMod.Ptype3,'val');
if Ptype1==1
    Ptype = 1;
elseif Ptype2==1
    Ptype = 2;
elseif Ptype3==1
    Ptype = 3;
end
if GCtype==1
    AS_ShowGCmatrix_residual(Ptype)
else
    AS_ShowGCmatrix_coefficient(Ptype)
end
uiwait(msgbox('Finish Load'));
end
function Return(varargin)
SelMod = varargin{3};
close(SelMod.fig);
ASBC_VIEW;
end

