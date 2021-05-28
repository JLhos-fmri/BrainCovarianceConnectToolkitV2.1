function BCCT_ShowMap_GUI
D.fig = figure('Name','Showing TFRZP maps',...            
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'color',[0.95 0.95 0.95],...
    'position',[0.2 0.3 0.6 0.4]);
movegui(D.fig,'center'); 
dires = which('BCCT_ShowMap_GUI.m');
[pth,nam,ext] = fileparts(dires);
D.pth = pth;

D.buttongroup1 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.45,0.8,0.4]);
D.buttongroup11 = uibuttongroup('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.01,2/3,0.98,1/3]);
D.TRFZ(1) = uicontrol('parent',D.buttongroup11,...
    'unit','norm',...
    'pos',[0.01,0,0.18,1],...
    'style','rad',...
    'string','T map');
D.TRFZ(2) = uicontrol('parent',D.buttongroup11,...
    'unit','norm',...
    'pos',[0.21,0,0.18,1],...
    'style','rad',...
    'string','F map');
D.TRFZ(3) = uicontrol('parent',D.buttongroup11,...
    'unit','norm',...
    'pos',[0.41,0,0.18,1],...
    'style','rad',...
    'string','R map');
D.TRFZ(4) = uicontrol('parent',D.buttongroup11,...
    'unit','norm',...
    'pos',[0.61,0,0.18,1],...
    'style','rad',...
    'string','Z map');
D.TRFZ(5) = uicontrol('parent',D.buttongroup11,...
    'unit','norm',...
    'pos',[0.81,0,0.18,1],...
    'style','rad',...
    'string','Other map');
D.TRFZmaptext = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.01,1/3,0.18,1/3],...
    'style','text',...
    'string',{'T/F/R/Z','Other maps'});
D.TFRZmapedit = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.2,1/3,0.39,1/3],...
    'style','edit',...
    'string','NULL');
D.TRFZmapsel = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.6,1/3,0.08,1/3],...
    'style','pushbutton',...
    'string','...');

D.TRFZmaptext1 = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.01,1/6,0.18,1/6],...
    'style','rad',...
    'string',{'with related pmap'});
D.TRFZmaptext2 = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.01,0,0.18,1/6],...
    'style','rad',...
    'string',{'without related pmap'});

D.TFRZmapedit_p = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.2,0,0.39,1/3],...
    'style','edit',...
    'string','NULL');
D.TRFZmapsel_p = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.6,0,0.08,1/3],...
    'style','pushbutton',...
    'string','...');
D.doftext = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.7,1/3,0.1,1/3],...
    'style','text',...
    'string',{'Degree of','Freedom'});
D.dofedit1 = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.8,1/3,0.1,1/3],...
    'style','edit',...
    'string','NULL');
D.dofedit2 = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.9,1/3,0.1,1/3],...
    'style','edit',...
    'string','NULL');

D.Int = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.7,0,0.07,1/3],...
    'style','pushbutton',...
    'string','Intensity');
D.IntE = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.77,0,0.07,1/3],...
    'style','text',...
    'string','1.7',...
    'horizontalalign','center');

D.Pva = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.85,0,0.07,1/3],...
    'style','pushbutton',...
    'string','Pval');

D.PvaE = uicontrol('parent',D.buttongroup1,...
    'unit','norm',...
    'pos',[0.92,0,0.07,1/3],...
    'style','text',...
    'string','0.05',...
    'horizontalalign','center');
%%
D.buttongroup2 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.325,0.8,0.1]);
D.Pmaptext = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.01,0.1,0.18,0.8],...
    'style','text',...
    'string','Pmap');
D.Pmapedit = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.2,0.1,0.39,0.8],...
    'style','edit',...
    'string','NULL');
D.Pmapsel = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.6,0.1,0.08,0.8],...
    'style','pushbutton',...
    'string','...');


D.PmapPval = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.7,0.1,0.14,0.8],...
    'style','pushbutton',...
    'string','Pval');

D.PmapPvalE = uicontrol('parent',D.buttongroup2,...
    'unit','norm',...
    'pos',[0.85,0.1,0.14,0.8],...
    'style','text',...
    'string','0.05',...
    'horizontalalign','center');
%%
D.buttongroup3 = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.2,0.8,0.1]);

D.BGtext = uicontrol('parent',D.buttongroup3,...
    'unit','norm',...
    'pos',[0.01,0.1,0.18,0.8],...
    'style','text',...
    'string','BackGround');
D.BGedit = uicontrol('parent',D.buttongroup3,...
    'unit','norm',...
    'pos',[0.2,0.1,0.69,0.8],...
    'style','edit',...
    'string','NULL');
D.BGsel = uicontrol('parent',D.buttongroup3,...
    'unit','norm',...
    'pos',[0.9,0.1,0.08,0.8],...
    'style','pushbutton',...
    'string','...');
set(D.BGsel,'callback',{@BGmapsel,D});
set(D.BGedit,'string',fullfile(pth,'mni_icbm152_t1.nii'));
%%
D.TRFsel = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.875,0.8,0.08]);

D.TRFselrad = uicontrol('parent',D.TRFsel,...
    'unit','norm',...
    'pos',[0.1,0.01,0.35,0.98],...
    'style','rad',...
    'string','TRFZmap');

D.Pselrad = uicontrol('parent',D.TRFsel,...
    'unit','norm',...
    'pos',[0.55,0.01,0.35,0.98],...
    'style','rad',...
    'string','Pmap');

set(D.TRFZ,'enable','on')
set(D.TRFZmaptext,'enable','on')
set(D.TFRZmapedit,'enable','on')
set(D.TRFZmapsel,'enable','on')
set(D.TRFZmaptext1,'enable','on');
set(D.TRFZmaptext2,'enable','on');
set(D.TFRZmapedit_p,'enable','on');
set(D.TRFZmapsel_p,'enable','on');
set(D.doftext,'enable','on');
set(D.dofedit1,'enable','on');
set(D.dofedit2,'enable','on');
set(D.Int,'enable','on');
set(D.IntE,'enable','on');
set(D.Pva,'enable','on');
set(D.PvaE,'enable','on');

set(D.Pmaptext,'enable','off')
set(D.Pmapedit,'enable','off')
set(D.Pmapsel,'enable','off')
set(D.PmapPval,'enable','off')
set(D.PmapPvalE,'enable','off')

set(D.TRFselrad,'callback',{@TRFsel,D})
set(D.Pselrad,'callback',{@Psel,D})

%%
D.showmod  = uibuttongroup('parent',D.fig,...
    'unit','norm',...
    'pos',[0.1,0.1,0.4,0.1]);
D.showboth = uicontrol('parent',D.showmod,...
    'unit','norm',...
    'pos',[0,0.1,1/3,0.8],...
    'style','rad',...
    'string','both');
D.showpos = uicontrol('parent',D.showmod,...
    'unit','norm',...
    'pos',[1/3,0.1,1/3,0.8],...
    'style','rad',...
    'string','pos only');
D.showneg = uicontrol('parent',D.showmod,...
    'unit','norm',...
    'pos',[2/3,0.1,1/3,0.8],...
    'style','rad',...
    'string','neg only');
%%
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
set(D.TRFZ(1),'callback',{@Tsel,D});
set(D.TRFZ(2),'callback',{@Fsel,D});
set(D.TRFZ(3),'callback',{@Rsel,D});
set(D.TRFZ(4),'callback',{@Zsel,D});
set(D.TRFZ(5),'callback',{@Othersel,D});
%%
set(D.TRFZmapsel,'callback',{@MapSel,D});
set(D.TRFZmapsel_p,'callback',{@TRF_PmapSel,D});
set(D.Pmapsel,'callback',{@PvaMapSel,D});
set(D.Int,'callback',{@IntSel,D});
set(D.Pva,'callback',{@PvaSel,D});
set(D.PmapPval,'callback',{@PmapPvalSel,D});

set(D.Exit,'callback',{@mapExit,D});
set(D.showbut,'callback',{@ShowButton,D});
end
function ShowButton(varargin)
D = varargin{3};
TRFZP = get(D.TRFselrad,'val');

if TRFZP
    Parameter.Mod = 1;
    TRFZ_Val = get(D.TRFZ,'val');
    if TRFZ_Val{1}; Parameter.Mod1.TRFZOther = 1;
    elseif TRFZ_Val{2}; Parameter.TRFZOther = 2;
    elseif TRFZ_Val{3}; Parameter.TRFZOther = 3;
    elseif TRFZ_Val{4}; Parameter.TRFZOther = 4;
    elseif TRFZ_Val{5}; Parameter.TRFZOther = 5;
    end
    TRFZmap = get(D.TFRZmapedit,'string');
    Parameter.Mod1.TRFZmap = TRFZmap;
    TRFZmap_p = get(D.TFRZmapedit_p,'string');
    Parameter.Mod1.TRFZmap_P = TRFZmap_p;
    Int = str2num(get(D.IntE,'string'));
    Parameter.Mod1.Int = Int;
    Pval = str2num(get(D.PvaE,'string'));
    Parameter.Mod1.Pval = Pval;
    PmapElab = get(D.TRFZmaptext1,'val');
    Parameter.Mod1.PmapElab = PmapElab;   
else
    Parameter.Mod = 2;    
    Pmap = get(D.Pmapedit,'string');
    Parameter.Mod2.Pmap = Pmap;
    Pval = str2num(get(D.PmapPvalE,'string'));
    Parameter.Mod2.Pval = Pval;    
end

bot = get(D.showboth,'val');
pos = get(D.showpos,'val');
neg = get(D.showneg,'val');
if bot
    Parameter.show = 0;
elseif pos
    Parameter.show = 1;
elseif neg
    Parameter.show = 2;
end
Parameter.pth = D.pth;
Parameter.BG = get(D.BGedit,'string');
if ~isempty(dir(fullfile(D.pth,'ShowParameter.mat')))
    delete(fullfile(D.pth,'ShowParameter.mat'));
end
save(fullfile(D.pth,'ShowParameter.mat'),'Parameter');
BCCT_showmap_subfun1(Parameter);
end


function BGmapsel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'background map selection');
PGseed = fullfile(PathName,FileName);
set(D.BGedit,'string',PGseed);
end

function mapExit(varargin)
D = varargin{3};
close(D.fig);
BCCT_VIEWmain;
end
function PmapPvalSel(varargin)
D = varargin{3};
Inputval = inputdlg('Please ENTER Pvals','Please ENTER Pvals',1,{'0.05'});
set(D.PmapPvalE,'string',Inputval{1});
end
function PvaSel(varargin)
D = varargin{3};
Inputval = inputdlg('Please ENTER Pvals','Please ENTER Pvals',1,{'0.05'});
set(D.PvaE,'string',Inputval{1});
pval = str2num(Inputval{1});
TRFZ_Val = get(D.TRFZ,'val');
if TRFZ_Val{1}     % T
    df1 = str2num(get(D.dofedit1,'string'));
    [Z,out] = AS_PtoTRF(pval,'T',df1,[]);
    set(D.IntE,'string',num2str(out));
elseif TRFZ_Val{2} % F
    df1 = str2num(get(D.dofedit1,'string'));
    df2 = str2num(get(D.dofedit2,'string'));
    [Z,out] = AS_PtoTRF(pval,'F',df1,df2);
    set(D.IntE,'string',num2str(out));
    
elseif TRFZ_Val{3} % R
    df1 = str2num(get(D.dofedit1,'string'));
    [Z,out] = AS_PtoTRF(pval,'R',df1,[]);
    set(D.IntE,'string',num2str(out));
    
elseif TRFZ_Val{4} % Z
    df1 = str2num(get(D.dofedit1,'string'));
    Zval = PtoZ(P);
    set(D.IntE,'string',num2str(Zval));
    
elseif TRFZ_Val{5} % Other
%     df1 = str2num(get(D.dofedit1,'string'));
%     pval = str2num(get(D.PvaE,'string'));
%     [Z,out] = AS_PtoTRF(pval,'T',df1,[]);
%     set(D.IntE,'string',num2str(out));
end
end
function IntSel(varargin)
D = varargin{3};

Inputval = inputdlg('Please ENTER Intensity values','Please ENTER Intensity values',1,{'1.7'});
set(D.IntE,'string',Inputval{1});
Intval = str2num(Inputval{1});
TRFZ_Val = get(D.TRFZ,'val');
if TRFZ_Val{1}     % T
    df1 = str2num(get(D.dofedit1,'string'));
    [Z,Pval] = AS_TFRtoZ(Intval,'T',df1,[]);
    set(D.PvaE,'string',num2str(Pval));
elseif TRFZ_Val{2} % F
    df1 = str2num(get(D.dofedit1,'string'));
    df2 = str2num(get(D.dofedit2,'string'));
    [Z,Pval] = AS_TFRtoZ(Intval,'F',df1,df2);
    set(D.PvaE,'string',num2str(Pval));
    
elseif TRFZ_Val{3} % R
    df1 = str2num(get(D.dofedit1,'string'));
    [Z,Pval] = AS_TFRtoZ(Intval,'R',df1,[]);
    set(D.PvaE,'string',num2str(Pval));
    
elseif TRFZ_Val{4} % Z
    Pval = normcdf(Intval);
    set(D.PvaE,'string',num2str(Pval));
    
elseif TRFZ_Val{5} % Other
%     df1 = str2num(get(D.dofedit1,'string'));
%     pval = str2num(get(D.PvaE,'string'));
%     [Z,out] = AS_PtoTRF(pval,'T',df1,[]);
%     set(D.IntE,'string',num2str(out));
end
end
function TRF_PmapSel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'P map selection');
PGseed = fullfile(PathName,FileName);
set(D.TFRZmapedit_p,'string',PGseed);
end
function PvaMapSel(varargin)
D = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'P map selection');
PGseed = fullfile(PathName,FileName);
set(D.Pmapedit,'string',PGseed);
end
function MapSel(varargin)
D = varargin{3};
TRFZ_Val = get(D.TRFZ,'val');
if TRFZ_Val{1}    
    [FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'T map selection');
    PGseed = fullfile(PathName,FileName);
    [v d] = Dynamic_read_dir_NIFTI(PGseed);
    desipinfo = v.descrip;
    inddes1 = find(desipinfo=='[');
    inddes2 = find(desipinfo==']');    
    if isempty(inddes1)
        dofans = inputdlg('dof of stat','dof of T stat',1,{'10'});
        dof = str2num(dofans{1});
    else
        dof = str2num(desipinfo(inddes1+1:inddes2-1));
    end
    set(D.dofedit1,'string',num2str(dof));
    [Z,out] = AS_PtoTRF(0.05,'T',dof,[]);
    out = round(out*1000)/1000;
    set(D.IntE,'string',num2str(out));
    set(D.TFRZmapedit,'string',PGseed);
elseif TRFZ_Val{2}
    [FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'F map selection');
    PGseed = fullfile(PathName,FileName);
    [v d] = Dynamic_read_dir_NIFTI(PGseed);
    desipinfo = v.descrip;
    inddes1 = find(desipinfo=='[');
    inddes2 = find(desipinfo==']');    
    if isempty(inddes1)
        dofans = inputdlg({'dof 1 of F stat','dof 2 of F stat'},'dof of F stat',1,{'2','10'});
        dof(1) = str2num(dofans{1});
        dof(2) = str2num(dofans{2});
    else
        dof = str2num(desipinfo(inddes1+1:inddes2-1));
    end
    set(D.dofedit1,'string',num2str(dof(1)));
    set(D.dofedit2,'string',num2str(dof(2)));
    [Z,out] = AS_PtoTRF(0.05,'F',dof(1),dof(2));
    out = round(out*1000)/1000;
    set(D.IntE,'string',num2str(out));
    set(D.TFRZmapedit,'string',PGseed);
elseif TRFZ_Val{3}
    [FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'R map selection');
    PGseed = fullfile(PathName,FileName);
    [v d] = Dynamic_read_dir_NIFTI(PGseed);
    desipinfo = v.descrip;
    inddes1 = find(desipinfo=='[');
    inddes2 = find(desipinfo==']');    
    if isempty(inddes1)
        dofans = inputdlg('dof of stat','dof of R stat',1,{'10'});
        dof = str2num(dofans{1});
    else
        dof = str2num(desipinfo(inddes1+1:inddes2-1));
    end
    set(D.dofedit1,'string',num2str(dof));
    [Z,out] = AS_PtoTRF(0.05,'R',dof,[]);
    out = round(out*1000)/1000;
    set(D.IntE,'string',num2str(out));
    set(D.TFRZmapedit,'string',PGseed);
elseif TRFZ_Val{4}
    [FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'Z map selection');
    PGseed = fullfile(PathName,FileName);
    Zval = AS_PtoZ(0.05);
    set(D.IntE,'string',num2str(Zval));
    set(D.TFRZmapedit,'string',PGseed);
elseif TRFZ_Val{5}
    [FileName,PathName,FilterIndex] = uigetfile({'*.nii';'*.img'},'Other map selection');
    PGseed = fullfile(PathName,FileName);
    set(D.TFRZmapedit,'string',PGseed);
end

end
function TRFsel(varargin)
D = varargin{3};
set(D.TRFselrad,'val',1);
set(D.Pselrad,'val',0);
% set(D.buttongroup1,'enable','on')
set(D.TRFZ,'enable','on')
set(D.TRFZmaptext,'enable','on')
set(D.TFRZmapedit,'enable','on')
set(D.TRFZmapsel,'enable','on')
set(D.TRFZmaptext1,'enable','on');
set(D.TRFZmaptext2,'enable','on');
set(D.TFRZmapedit_p,'enable','on');
set(D.TRFZmapsel_p,'enable','on');
set(D.doftext,'enable','on');
set(D.dofedit1,'enable','on');
set(D.dofedit2,'enable','on');
set(D.Int,'enable','on');
set(D.IntE,'enable','on');
set(D.Pva,'enable','on');
set(D.PvaE,'enable','on');


% set(D.buttongroup2,'enable','off')
set(D.Pmaptext,'enable','off')
set(D.Pmapedit,'enable','off')
set(D.Pmapsel,'enable','off')
set(D.PmapPval,'enable','off')
set(D.PmapPvalE,'enable','off')
end
function Psel(varargin)
D = varargin{3};
set(D.TRFselrad,'val',0);
set(D.Pselrad,'val',1);
% set(D.buttongroup1,'enable','on')
set(D.TRFZ,'enable','off')
set(D.TRFZmaptext,'enable','off')
set(D.TFRZmapedit,'enable','off')
set(D.TRFZmapsel,'enable','off')
set(D.TRFZmaptext1,'enable','off');
set(D.TRFZmaptext2,'enable','off');
set(D.TFRZmapedit_p,'enable','off');
set(D.TRFZmapsel_p,'enable','off');
set(D.doftext,'enable','off');
set(D.dofedit1,'enable','off');
set(D.dofedit2,'enable','off');
set(D.Int,'enable','off');
set(D.IntE,'enable','off');
set(D.Pva,'enable','off');
set(D.PvaE,'enable','off');

% set(D.buttongroup2,'enable','off')
set(D.Pmaptext,'enable','on')
set(D.Pmapedit,'enable','on')
set(D.Pmapsel,'enable','on')
set(D.PmapPval,'enable','on')
set(D.PmapPvalE,'enable','on')
end
function Tsel(varargin)
D = varargin{3};
set(D.dofedit1,'enable','on')
set(D.dofedit2,'enable','off')
set(D.Pva,'enable','on');
set(D.PvaE,'enable','on');
end
function Fsel(varargin)
D = varargin{3};
set(D.dofedit1,'enable','on')
set(D.dofedit2,'enable','on')
set(D.Pva,'enable','on');
set(D.PvaE,'enable','on');
end
function Rsel(varargin)
D = varargin{3};
set(D.dofedit1,'enable','on')
set(D.dofedit2,'enable','off')
set(D.Pva,'enable','on');
set(D.PvaE,'enable','on');
end
function Zsel(varargin)
D = varargin{3};
set(D.dofedit1,'enable','off')
set(D.dofedit2,'enable','off')
set(D.Pva,'enable','on');
set(D.PvaE,'enable','on');
end
function Othersel(varargin)
D = varargin{3};
set(D.dofedit1,'enable','off')
set(D.dofedit2,'enable','off')
set(D.Pva,'enable','off');
set(D.PvaE,'enable','off');
end
