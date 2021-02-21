function BCCT_ShowFCmatrixGUI_GCA
Hsize = get(0,'screensize');
MIDPOINT = [Hsize(3)/2,Hsize(4)/2];
Asize = [100*3,100+40];
MaxSIZE = [Hsize(3) Hsize(4)]*0.4;
factor = MaxSIZE./Asize;
factornew = min(factor);
POSSIZE = Asize*factornew;
Hshow.fig = figure('position',[MIDPOINT(1)-POSSIZE(1)/2,MIDPOINT(2)-POSSIZE(2)/2,POSSIZE(1),POSSIZE(2)],'name','CaSCN Matrix Show');

Hshow.IO = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.1 0.7 0.8 0.2]);
Hshow.IO_inputpb = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.05 0.6 0.1 0.3],'style','text','string','Inputdir');
Hshow.IO_inputsel = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.85 0.6 0.1 0.3],'style','pushbutton','string','Select');
Hshow.IO_inputed = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.2 0.6 0.6 0.3],'style','edit');
Hshow.IO_Outputpb = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.05 0.1 0.1 0.3],'style','text','string','Outputdir');
Hshow.IO_Outputsel = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.85 0.1 0.1 0.3],'style','pushbutton','string','Select');
pwdpath = pwd;
Hshow.IO_Outputed = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.2 0.1 0.6 0.3],'style','edit','string',pwdpath);

Hshow.COLORMAP = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.1 0.55 0.8 0.1]);
Hshow.Colmap1 = uicontrol('parent',Hshow.COLORMAP,'units','norm','pos',[0.1 0.05 0.2 0.9],'style','rad','string','Colormap:AFNI');
Hshow.Colmap2 = uicontrol('parent',Hshow.COLORMAP,'units','norm','pos',[0.4 0.05 0.2 0.9],'style','rad','string','Colormap:RdBu');
Hshow.Colmap3 = uicontrol('parent',Hshow.COLORMAP,'units','norm','pos',[0.7 0.05 0.2 0.9],'style','rad','string','Colormap:RdYlBu');
set(Hshow.Colmap1,'val',1);
set(Hshow.Colmap2,'val',0);
set(Hshow.Colmap3,'val',0);

Hshow.DOFg = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.1 0.4 0.2 0.1]);
Hshow.res = uicontrol('parent',Hshow.DOFg,'unit','norm','pos',[0.1 0.1 0.35 0.8],'style','rad','string','res')
Hshow.coef = uicontrol('parent',Hshow.DOFg,'unit','norm','pos',[0.55 0.1 0.35 0.8],'style','rad','string','coef')
% Hshow.DOF = uicontrol('parent',Hshow.DOFg,'units','norm','pos',[0.1 0.05 0.4 0.9],'style','text','string','Degree of Freedom');
% Hshow.DOFE = uicontrol('parent',Hshow.DOFg,'units','norm','pos',[0.5 0.05 0.4 0.9],'style','edit');

Hshow.TYPE = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.4 0.4 0.2 0.1]);
Hshow.Type1 = uicontrol('parent',Hshow.TYPE,'units','norm','pos',[0.1 0.05 0.2 0.9],'style','rad','string','Parameter');
Hshow.Type2 = uicontrol('parent',Hshow.TYPE,'units','norm','pos',[0.5 0.05 0.4 0.9],'style','rad','string','Perm P');

Hshow.PVAL = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.7 0.4 0.2 0.1]);
Hshow.Pval = uicontrol('parent',Hshow.PVAL,'units','norm','pos',[0.1 0.05 0.4 0.9],'style','text','string','P value');
Hshow.PvalE = uicontrol('parent',Hshow.PVAL,'units','norm','pos',[0.5 0.05 0.4 0.9],'style','edit');

Hshow.PTYPE = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.1 0.25 0.8 0.1]);
Hshow.PTYPE1 = uicontrol('parent',Hshow.PTYPE,'units','norm','pos',[0.075 0.05 0.2 0.9],'style','rad','string','Uncorrected');
Hshow.PTYPE2 = uicontrol('parent',Hshow.PTYPE,'units','norm','pos',[0.3 0.05 0.2 0.9],'style','rad','string','FPA');
Hshow.PTYPE3 = uicontrol('parent',Hshow.PTYPE,'units','norm','pos',[0.525 0.05 0.2 0.9],'style','rad','string','FDR');
Hshow.PTYPE4 = uicontrol('parent',Hshow.PTYPE,'units','norm','pos',[0.75 0.05 0.2 0.9],'style','rad','string','Bonf');

Hshow.ORDERS = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.1 0.1 0.2 0.1]);
Hshow.order = uicontrol('parent',Hshow.ORDERS,'units','norm','pos',[0.1 0.05 0.3 0.9],'style','text','string','Orders');
Hshow.orderE = uicontrol('parent',Hshow.ORDERS,'units','norm','pos',[0.4 0.05 0.5 0.9],'style','edit');
Hshow.UNDERLINE = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.4 0.1 0.2 0.1]);
Hshow.underline = uicontrol('parent',Hshow.UNDERLINE,'units','norm','pos',[0.1 0.05 0.3 0.9],'style','text','string','Underline index');
Hshow.underlineE = uicontrol('parent',Hshow.UNDERLINE,'units','norm','pos',[0.4 0.05 0.5 0.9],'style','edit');

Hshow.EVAL = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.7 0.1 0.2 0.1]);
Hshow.Show = uicontrol('parent',Hshow.EVAL,'units','norm','pos',[0.1 0.05 0.35 0.9],'style','pushbutton','string','Show');
Hshow.Exit = uicontrol('parent',Hshow.EVAL,'units','norm','pos',[0.55 0.05 0.35 0.9],'style','pushbutton','string','Exit');
set(Hshow.IO_inputed,'callback',{@HSMATRIXinputedit,Hshow});
set(Hshow.IO_inputsel,'callback',{@HSMATRIXinput,Hshow});
set(Hshow.IO_Outputsel,'callback',{@HSMATRIXoutput,Hshow});
set(Hshow.Show,'callback',{@HSMATRIXSHOW,Hshow});
set(Hshow.Exit,'callback',{@HSMATRIXEXIT,Hshow});

end

function HSMATRIXinputedit(varargin)
Hshow = varargin{3};
patnam = get(Hshow.IO_inputed,'string');
tempinfo = load(patnam);
labval = get(Hshow.res,'val');
GCAres = tempinfo.GCAres;
if labval % res
    Outshow.GCAlab = 1;
    Outshow.res.Pval_x2y = GCAres.GCA_res_x2y_pval;
    Outshow.res.Pval_y2x = GCAres.GCA_res_y2x_pval;
    Outshow.res.Orig_x2y = GCAres.GCA_res_x2y;
    Outshow.res.Orig_y2x = GCAres.GCA_res_y2x;
    Outshow.res.Trans_x2y = GCAres.GCA_res_x2y_trans;
    Outshow.res.Trans_y2x = GCAres.GCA_res_y2x_trans;
    if isfield(GCAres,'GCA_res_x2y_permpval')
        Outshow.permlab = 1;
        Outshow.res.Permp_Orig_x2y = GCAres.GCA_res_x2y_permpval;
        Outshow.res.Permp_Orig_y2x = GCAres.GCA_res_y2x_permpval;
        Outshow.res.Permp_Trans_x2y = GCAres.GCA_res_x2y_trans_permpval;
        Outshow.res.Permp_Trans_y2x = GCAres.GCA_res_y2x_trans_permpval;
        Outshow.res.Permp_Fx2y = GCAres.GCA_res_fx2y_permpval;
        Outshow.res.Fx2y = GCAres.GCA_res_fx2y;
    end
else % coef
    Outshow.GCAlab = 2;
    
end

end
function HSMATRIXinput(varargin)
Hshow = varargin{3};

[nam,path,ext] = uigetfile('*.mat','Matrix of Connection');
patnam = fullfile(path,nam);
labval = get(Hshow.res,'val');
set(Hshow.IO_inputed,'string',patnam);

if labval % res
    
else % coef
    
end
[pat,nam0,ext0] = fileparts(which('BCCT_ShowFCmatrixGUI.m'));
if isempty(dir([pat,filesep,'MatrixOutshow',filesep,'SetUpinfo.mat']))
    save([pat,filesep,'MatrixOutshow',filesep,'SetUpinfo.mat'],'ShowMatinfo');
else
    delete([pat,filesep,'MatrixOutshow',filesep,'SetUpinfo.mat']);
    save([pat,filesep,'MatrixOutshow',filesep,'SetUpinfo.mat'],'ShowMatinfo');
end

end

function HIOSELTYP1(varargin)
Hio = varargin{3};
set(Hio.Type1,'val',1);
set(Hio.Type2,'val',0);
end
function HIOSELTYP2(varargin)
Hio = varargin{3};
set(Hio.Type1,'val',0);
set(Hio.Type2,'val',1);
end
function HIOOK(varargin)
Hio = varargin{3};
Types = get(Hio.Type1,'val');
DOFV = get(Hio.DFE,'string');
NROIS = get(Hio.NROIE,'string');
nrois = str2num(NROIS);
matnam = get(Hio.MatnameE,'string');
[pat,nam,ext] = fileparts(which('BCCT_ShowFCmatrixGUI.m'));
if isempty(dir([pat,filesep,'MatrixOutshow']))
    mkdir([pat,filesep,'MatrixOutshow']);
else
    rmdir([pat,filesep,'MatrixOutshow'],'s');
    mkdir([pat,filesep,'MatrixOutshow']);
end
INFOS.matnam = matnam;
INFOS.Types = Types;
INFOS.DOFV = DOFV;
INFOS.nrois = nrois;
save(fullfile(pat,'MatrixOutshow','INFOS.mat'),'INFOS');
close(Hio.fig);
% uiwait(msgbox('Info Load OK!'));
end

function HSMATRIXoutput(varargin)
Hshow = varargin{3};
Path = uigetdir(pwd,'Output for pic and other information');
set(Hshow.IO_Outputed,'string',Path);
end
function HSMATRIXSHOW(varargin)
Hshow = varargin{3};
Outputdir = get(Hshow.IO_Outputed,'string');
Inputdir = get(Hshow.IO_inputed,'string');

[pat,nam0,ext0] = fileparts(which('BCCT_ShowFCmatrixGUI.m'));
load([pat,filesep,'MatrixOutshow',filesep,'SetUpinfo.mat'])
Types = ShowMatinfo.Type;
DF_E = ShowMatinfo.DF_E;
% save temps
Colmaptype1 = get(Hshow.Colmap1,'val');
Colmaptype2 = get(Hshow.Colmap2,'val');
Colmaptype3 = get(Hshow.Colmap3,'val');
if Colmaptype1
    Coltype = 1;
elseif Colmaptype2
    Coltype = 2;
else
    Coltype = 3;
end
Datatype1 = get(Hshow.Type1,'val');
Datatype2 = get(Hshow.Type2,'val');
if Datatype1
    TYPES = 1;
else
    TYPES = 2;
end
PVALT1 = get(Hshow.PTYPE1,'val');
PVALT2 = get(Hshow.PTYPE2,'val');
PVALT3 = get(Hshow.PTYPE3,'val');
PVALT4 = get(Hshow.PTYPE4,'val');
if PVALT1 
    PVALT = 1;
elseif PVALT2 
    PVALT = 2;
elseif PVALT3 
    PVALT = 3;
else
    PVALT = 4;
end
Pvale = get(Hshow.PvalE,'string');
Pvalused = str2num(Pvale);
enhanindedit = get(Hshow.underlineE,'string');
enhanind = str2num(enhanindedit);
ORDERSHOW = get(Hshow.orderE,'string');
orders = str2num(ORDERSHOW);
if Types==1 % FCmat
    if TYPES~=Types
        error('wrong input data type, it is R matrix of ASBC');
    end    
    AS_MatrixFCshow_bat_GUI(Inputdir,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders);
elseif Types==2 % Pmap    
    if TYPES~=Types
        error('wrong input data type, it is Perm P matrix of ASBC');
    end
    AS_MatrixShowFinalResPermVer_bat_GUI(Inputdir,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders)
elseif Types==31
    R = ShowMatinfo.R;
%     DOF = str2num(DF_E);
    AS_MatrixFCshow_other_GUI(R,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders,DF_E)
elseif Types==32
    Pval = ShowMatinfo.Pval;
    AS_MatrixShowFinalResPermVer_other_GUI(Pval,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders)
end
end
function HSMATRIXEXIT(varargin)
Hshow = varargin{3};
close(Hshow.fig);
ASBC_VIEW;
end