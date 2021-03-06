function BCCT_ShowGCMap_GUI_res
D.fig = figure('Name','Showing Res-based GCA maps',...            
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'color',[0.95 0.95 0.95],...
    'position',[0.2 0.3 0.6 0.4]);
movegui(D.fig,'center'); 
dires = which('BCCT_ShowGCMap_GUI_res.m');
[pth,nam,ext] = fileparts(dires);
D.pth = pth;


D.buttongroup1 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.8,0.8,0.1]);
D.Maptype(1) = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.1,0.1,0.35,0.8],...
    'style','rad',...
    'string','Res map');
D.Maptype(2) = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.55,0.1,0.35,0.8],...
    'style','rad',...
    'string','Perm Pmap');
%%
D.buttongroup2 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.475,0.8,0.3]);
D.coefmap_text = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.1,0.65,0.1,0.3],...
    'style','text',...
    'string',{'ResMap','(Res/Trans)'});

D.coefmap_edit = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.2,0.65,0.5,0.3],...
    'style','edit',...
    'string','null');

D.coefmap_sel = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.7,0.65,0.1,0.3],...
    'style','pushbutton',...
    'string','...');

D.coefmap_p_text = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.1,0.35,0.1,0.3],...
    'style','rad',...
    'string','Pmap');

D.coefmap_p_edit = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.2,0.35,0.5,0.3],...
    'style','edit',...
    'string','null');

D.coefmap_p_sel = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.7,0.35,0.1,0.3],...
    'style','pushbutton',...
    'string','...');

D.coefmap_perm_text = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.1,0.05,0.1,0.3],...
    'style','rad',...
    'string','Perm map',...
    'enable','off');

D.coefmap_perm_edit = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.2,0.05,0.5,0.3],...
    'style','edit',...
    'string','null',...
    'enable','off');

D.coefmap_perm_sel = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.7,0.05,0.1,0.3],...
    'style','pushbutton',...
    'string','...',...
    'enable','off');

D.coefmap_thr_pb = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.85,0.6,0.1,0.3],...
    'style','pushbutton',...
    'string','Pvalue',...
    'enable','off');
D.coefmap_thr_text = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.85,0.1,0.1,0.3],...
    'style','text',...
    'string','0.05',...
    'enable','off');
set(D.coefmap_sel,'callback',{@selcoefmap,D});
set(D.coefmap_p_sel,'callback',{@selcoefmap_p,D});
set(D.coefmap_perm_sel,'callback',{@selcoefmap_perm,D});
set(D.coefmap_p_text,'callback',{@sel_pval1,D});
set(D.coefmap_perm_text,'callback',{@sel_pval2,D});
set(D.coefmap_thr_pb,'callback',{@coef_pthr,D})

%%
D.buttongroup3 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.35,0.8,0.1]);

D.perm_p_text = uicontrol('parent',D.buttongroup3,...
    'unit','norm',...
    'pos',[0.1,0.1,0.1,0.8],...
    'style','text',...
    'string','Perm map',...
    'enable','off');

D.perm_p_edit = uicontrol('parent',D.buttongroup3,...
    'unit','norm',...
    'pos',[0.2,0.1,0.5,0.8],...
    'style','edit',...
    'string','null',...
    'enable','off');

D.perm_p_sel = uicontrol('parent',D.buttongroup3,...
    'unit','norm',...
    'pos',[0.7,0.1,0.1,0.8],...
    'style','pushbutton',...
    'string','...',...
    'enable','off');
D.perm_p_thr_pb = uicontrol('parent',D.buttongroup3,...
    'unit','norm',...
    'pos',[0.8,0.1,0.1,0.8],...
    'style','pushbutton',...
    'string','Pvalue',...
    'enable','off');
D.perm_p_thr_text = uicontrol('parent',D.buttongroup3,...
    'unit','norm',...
    'pos',[0.9,0.1,0.1,0.8],...
    'style','text',...
    'string','0.05',...
    'enable','off');
set(D.perm_p_sel,'callback',{@perm_psel,D});
set(D.perm_p_thr_pb,'callback',{@perm_pthr,D});
%%
D.buttongroup4 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.225,0.8,0.1]);

D.BG_text = uicontrol('parent',D.buttongroup4,...
    'unit','norm',...
    'pos',[0.1,0.1,0.1,0.8],...
    'style','text',...
    'string',{'BackGround', 'map'});

D.BG_edit = uicontrol('parent',D.buttongroup4,...
    'unit','norm',...
    'pos',[0.2,0.1,0.6,0.8],...
    'style','edit',...
    'string','null');

D.BG_sel = uicontrol('parent',D.buttongroup4,...
    'unit','norm',...
    'pos',[0.8,0.1,0.1,0.8],...
    'style','pushbutton',...
    'string','...');
set(D.BG_edit,'string',fullfile(pth,'mni_icbm152_t1.nii'));
set(D.BG_sel,'callback',{@BGselfun,D});

%%
D.buttongroup5 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.1,0.4,0.1]);
D.showboth = uicontrol('parent',D.buttongroup5,...
    'unit','norm',...
    'pos',[0,0.1,1/3,0.8],...
    'style','rad',...
    'string','both');
D.showpos = uicontrol('parent',D.buttongroup5,...
    'unit','norm',...
    'pos',[1/3,0.1,1/3,0.8],...
    'style','rad',...
    'string','pos only');
D.showneg = uicontrol('parent',D.buttongroup5,...
    'unit','norm',...
    'pos',[2/3,0.1,1/3,0.8],...
    'style','rad',...
    'string','neg only');

D.showbut = uicontrol('parent',D.fig,...
    'unit','norm',...
    'pos',[0.6,0.1,0.15,0.1],...
    'style','pushbutton',...
    'string','Show');
D.Exit = uicontrol('parent',D.fig,...
    'unit','norm',...
    'pos',[0.75,0.1,0.15,0.1],...
    'style','pushbutton',...
    'string','Exit');


set(D.Maptype,'callback',{@selMapmod,D});
set(D.showbut,'callback',{@calcu,D});
set(D.Exit,'callback',{@Exitf,D});
end

function calcu(varargin)
D = varargin{3};
CPlab = get(D.Maptype(1),'val'); % 1 for coef, 2 for perm p
PNlab(1) = get(D.showboth,'val');
PNlab(2) = get(D.showpos,'val');
PNlab(3) = get(D.showneg,'val');
coefptype = get(D.coefmap_p_text,'val');
coefmap = get(D.coefmap_edit,'string');
coefpmap = get(D.coefmap_p_edit,'string');
coefpermpmap = get(D.coefmap_perm_edit,'string');
permpmap = get(D.perm_p_edit,'string');
coefpthr = get(D.coefmap_thr_text,'string');
permpthr = get(D.perm_p_thr_text,'string');
bgmap = get(D.BG_edit,'string');

Parameter.CPlab = CPlab;
Parameter.PNlab = PNlab;
Parameter.coefmap = coefmap;
Parameter.coefpmap = coefpmap;
Parameter.coefpermpmap = coefpermpmap;
Parameter.coefpthr = str2num(coefpthr);
Parameter.permpthr = str2num(permpthr);
Parameter.bgmap = bgmap;
Parameter.permpmap = permpmap;
Parameter.coefptype = coefptype;

save(fullfile(D.pth,'ShowParameter.mat'),'Parameter');
BCCT_showgcmap_subfun2(Parameter)
end
function Exitf(varargin)
D = varargin{3};
close(D.fig);
BCCT_VIEWmain;
end
function BGselfun(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'BackGround');
PGseed = fullfile(PathName,FileName);
set(D.BG_edit,'string',PGseed);
end
function perm_psel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'Permutation Pmap(mainly for group)');
PGseed = fullfile(PathName,FileName);
set(D.coefmap_edit,'string',PGseed);
end
function perm_pthr(varargin)
D = varargin{3};
Inputval = inputdlg('Please ENTER Pvals','Please ENTER Pvals',1,{'0.05'});
set(D.perm_p_thr_text,'string',Inputval{1});
end
function coef_pthr(varargin)
D = varargin{3};
Inputval = inputdlg('Please ENTER Pvals','Please ENTER Pvals',1,{'0.05'});
set(D.coefmap_thr_text,'string',Inputval{1});
end
function selcoefmap(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'coefmap (coef/T/Z)');
PGseed = fullfile(PathName,FileName);
set(D.coefmap_edit,'string',PGseed);
end
function selcoefmap_p(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'coefmap (pmap)');
PGseed = fullfile(PathName,FileName);
set(D.coefmap_p_edit,'string',PGseed);
end
function selcoefmap_perm(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'coefmap (perm p)');
PGseed = fullfile(PathName,FileName);
set(D.coefmap_perm_edit,'string',PGseed);
end
function selMapmod(varargin)
D = varargin{3};
val = get(D.Maptype(1),'val');
if val
   set(D.perm_p_edit,'enable','off');
   set(D.perm_p_sel,'enable','off');
   set(D.perm_p_text,'enable','off'); 
   set(D.perm_p_thr_pb,'enable','off');
   set(D.perm_p_thr_text,'enable','off');
   set(D.coefmap_edit,'enable','on');
   set(D.coefmap_sel,'enable','on');
   set(D.coefmap_text,'enable','on');
   set(D.coefmap_p_edit,'enable','on');
   set(D.coefmap_p_sel,'enable','on');
   set(D.coefmap_p_text,'enable','on');
   set(D.coefmap_perm_edit,'enable','on');
   set(D.coefmap_perm_sel,'enable','on');
   set(D.coefmap_perm_text,'enable','on');
   set(D.coefmap_thr_pb,'enable','on');
   set(D.coefmap_thr_text,'enable','on');   
else
   set(D.perm_p_edit,'enable','on');
   set(D.perm_p_sel,'enable','on');
   set(D.perm_p_text,'enable','on'); 
   set(D.perm_p_thr_pb,'enable','on');
   set(D.perm_p_thr_text,'enable','on');
   set(D.coefmap_edit,'enable','off');
   set(D.coefmap_sel,'enable','off');
   set(D.coefmap_text,'enable','off');
   set(D.coefmap_p_edit,'enable','off');
   set(D.coefmap_p_sel,'enable','off');
   set(D.coefmap_p_text,'enable','off');
   set(D.coefmap_perm_edit,'enable','off');
   set(D.coefmap_perm_sel,'enable','off');
   set(D.coefmap_perm_text,'enable','off');
   set(D.coefmap_thr_pb,'enable','off');
   set(D.coefmap_thr_text,'enable','off');    
end
end
function sel_pval1(varargin)
D = varargin{3};
set(D.coefmap_p_edit,'enable','on');
set(D.coefmap_p_sel,'enable','on');
set(D.coefmap_perm_edit,'enable','off');
set(D.coefmap_perm_sel,'enable','off');
end
function sel_pval2(varargin)
D = varargin{3};
set(D.coefmap_p_edit,'enable','off');
set(D.coefmap_p_sel,'enable','off');
set(D.coefmap_perm_edit,'enable','on');
set(D.coefmap_perm_sel,'enable','on');
end

