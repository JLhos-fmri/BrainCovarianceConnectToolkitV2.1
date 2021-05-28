function BCCT_ShowGCmatrix_GUI
Hsize = get(0,'screensize');
MIDPOINT = [Hsize(3)/2,Hsize(4)/2];
Asize = [100*3,100+40];
MaxSIZE = [Hsize(3) Hsize(4)]*0.4;
factor = MaxSIZE./Asize;
factornew = min(factor);
POSSIZE = Asize*factornew;
Hshow.fig = figure('position',[MIDPOINT(1)-POSSIZE(1)/2,MIDPOINT(2)-POSSIZE(2)/2,POSSIZE(1),POSSIZE(2)],'name','Connectivity Matrix Show');

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

% Hshow.DOFg = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.1 0.4 0.2 0.1]);
% Hshow.DOF = uicontrol('parent',Hshow.DOFg,'units','norm','pos',[0.1 0.05 0.4 0.9],'style','text','string','Degree of Freedom');
% Hshow.DOFE = uicontrol('parent',Hshow.DOFg,'units','norm','pos',[0.5 0.05 0.4 0.9],'style','edit');

Hshow.TYPE = uibuttongroup('parent',Hshow.fig,'units','norm','pos',[0.1 0.4 0.5 0.1]);
Hshow.Type1 = uicontrol('parent',Hshow.TYPE,'units','norm','pos',[0.05 0.05 0.3 0.9],'style','rad','string','Res');
Hshow.Type3 = uicontrol('parent',Hshow.TYPE,'units','norm','pos',[0.35 0.05 0.3 0.9],'style','rad','string','Coef');
Hshow.Type2 = uicontrol('parent',Hshow.TYPE,'units','norm','pos',[0.65 0.05 0.3 0.9],'style','rad','string','Perm P');

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

DatType(1,1) = get(Hshow.Type1,'val'); % res
DatType(2,1) = get(Hshow.Type3,'val'); % coef
DatType(3,1) = get(Hshow.Type2,'val'); % PermP

if DatType(2,1) % coef
    ResGCA = tempinfo.GCAres;
    Ord = size(ResGCA,2);
    if ~isfield(ResGCA(1),'GCA_coef_x2y')
        error('Not CaSCN results of Coef, please re-select');
    end
    ShowMatinfo.coef = ResGCA;
    ShowMatinfo.coeford = Ord;
    nrois = size(ResGCA(1).GCA_coef_x2y,1);
    ShowMatinfo.Type = 2;
elseif DatType(1,1) % res
    ResGCA = tempinfo.GCAres;
    if ~isfield(ResGCA(1),'GCA_res_x2y')
        error('Not CaSCN results of res, please re-select');
    end
    ShowMatinfo.res = ResGCA; 
    nrois = size(ResGCA(1).GCA_res_x2y,1);
    ShowMatinfo.Type = 1;
else % Perm Group Comp
    ResGCA = tempinfo.GroupCompPerm;
    if isfield(ResGCA(1),'res')
        ShowMatinfo.G2resperm = ResGCA.res;
        nrois = size(ResGCA.res.P_mat1,1);
    ShowMatinfo.Type = 3;
    elseif isfield(ResGCA(1),'coef')
        ShowMatinfo.G2coefperm = ResGCA.coef;
        nrois = size(ResGCA.coef.P_mat1,1);
    ShowMatinfo.Type = 4;
    else
        error('No group comparision results of CaSCN, please re-select');
    end
end
set(Hshow.orderE,'string',num2str(1:nrois));
set(Hshow.underlineE,'string',num2str([1,nrois]));


[pat,nam0,ext0] = fileparts(which('BCCT_ShowGCmatrixGUI.m'));
if ~isdir([pat,filesep,'MatrixOutshowGC',filesep])
    mkdir([pat,filesep,'MatrixOutshowGC',filesep]);
end
if isempty(dir([pat,filesep,'MatrixOutshowGC',filesep,'SetUpinfo.mat']))
    save([pat,filesep,'MatrixOutshowGC',filesep,'SetUpinfo.mat'],'ShowMatinfo');
else
    delete([pat,filesep,'MatrixOutshowGC',filesep,'SetUpinfo.mat']);
    save([pat,filesep,'MatrixOutshowGC',filesep,'SetUpinfo.mat'],'ShowMatinfo');
end
end
function HSMATRIXinput(varargin)
Hshow = varargin{3};
[nam,path,ext] = uigetfile('*.mat','Matrix of Connection');
patnam = fullfile(path,nam);
set(Hshow.IO_inputed,'string',patnam);
tempinfo = load(patnam);
DatType(1,1) = get(Hshow.Type1,'val'); % res
DatType(2,1) = get(Hshow.Type3,'val'); % coef
DatType(3,1) = get(Hshow.Type2,'val'); % PermP

if DatType(1,1) % res
    ResGCA = tempinfo.GCAres;
    if ~isfield(ResGCA(1),'GCA_res_x2y')
        error('Not CaSCN results of res, please re-select');
    end
    ShowMatinfo.res = ResGCA; 
    nrois = size(ResGCA(1).GCA_res_x2y,1);
    ShowMatinfo.Type = 1;
elseif DatType(2,1) % coef
    ResGCA = tempinfo.GCAres;
    Ord = size(ResGCA,2);
    if ~isfield(ResGCA(1),'GCA_coef_x2y')
        error('Not CaSCN results of Coef, please re-select');
    end
    ShowMatinfo.coef = ResGCA;
    ShowMatinfo.coeford = Ord;
    nrois = size(ResGCA(1).GCA_coef_x2y,1);
    ShowMatinfo.Type = 2;
else % Perm Group Comp
    ResGCA = tempinfo.GroupCompPerm;
    if isfield(ResGCA(1),'res')
        ShowMatinfo.G2resperm = ResGCA.res;
        nrois = size(ResGCA.res.P_mat1,1);
        ShowMatinfo.Type = 3;
    elseif isfield(ResGCA(1),'coef')
        ShowMatinfo.G2coefperm = ResGCA.coef;
        nrois = size(ResGCA.coef.P_mat1,1);
        ShowMatinfo.Type = 4;
    else
        error('No group comparision results of CaSCN, please re-select');
    end
end
set(Hshow.orderE,'string',num2str(1:nrois));
set(Hshow.underlineE,'string',num2str([1,nrois]));

[pat,nam0,ext0] = fileparts(which('BCCT_ShowGCmatrix_GUI.m'));

if ~isdir([pat,filesep,'MatrixOutshowGC',filesep])
    mkdir([pat,filesep,'MatrixOutshowGC',filesep]);
end
if isempty(dir([pat,filesep,'MatrixOutshowGC',filesep,'SetUpinfo.mat']))
    save([pat,filesep,'MatrixOutshowGC',filesep,'SetUpinfo.mat'],'ShowMatinfo');
else
    delete([pat,filesep,'MatrixOutshowGC',filesep,'SetUpinfo.mat']);
    save([pat,filesep,'MatrixOutshowGC',filesep,'SetUpinfo.mat'],'ShowMatinfo');
end

end
%%
function HSMATRIXoutput(varargin)
Hshow = varargin{3};
Path = uigetdir(pwd,'Output for pic and other information');
set(Hshow.IO_Outputed,'string',Path);
end
function HSMATRIXSHOW(varargin)
Hshow = varargin{3};
Outputdir = get(Hshow.IO_Outputed,'string');
Inputdir = get(Hshow.IO_inputed,'string');

[pat,nam0,ext0] = fileparts(which('BCCT_ShowGCmatrix_GUI.m'));
load([pat,filesep,'MatrixOutshowGC',filesep,'SetUpinfo.mat']);
Types = ShowMatinfo.Type;
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
if Datatype1 %res
    TYPES = 1;
elseif Datatype2 % perm
    TYPES = 3;
else % coef
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
if Types==1 % res
    if TYPES~=Types
        error('wrong input data type, it is residual matrix of CaSCN');
    end    
    RESshow = ShowMatinfo.res;
    BCCT_MatrixGCshow_res_GUI(RESshow,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders)
elseif Types==2 % coef
    if TYPES~=Types
        error('wrong input data type, it is coeficience matrix of CaSCN');
    end
    RESshow = ShowMatinfo.coef;
    BCCT_MatrixGCshow_coef_GUI(RESshow,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders)
    
%     AS_MatrixShowFinalResPermVer_bat_GUI(Inputdir,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders)
elseif Types==3 % res perm 
    Pvalu = ShowMatinfo.G2resperm;
%     save test1
%     DOF = str2num(DF_E);
    Pval = Pvalu.P_mat1;
    BCCT_MatrixShowFinalResPermVer_other_GUI(Pval,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders,'GroupComp-CaSCN-res-X2Y')
    Pval = Pvalu.P_mat2;
    BCCT_MatrixShowFinalResPermVer_other_GUI(Pval,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders,'GroupComp-CaSCN-res-Y2X')
    Pval = Pvalu.P_mat3;
    BCCT_MatrixShowFinalResPermVer_other_GUI(Pval,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders,'GroupComp-CaSCN-res-X2Y-TRANS')
    Pval = Pvalu.P_mat4;
    BCCT_MatrixShowFinalResPermVer_other_GUI(Pval,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders,'GroupComp-CaSCN-res-Y2X-TRANS')
    Pval = Pvalu.P_mat5;
    BCCT_MatrixShowFinalResPermVer_other_GUI(Pval,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders,'GroupComp-CaSCN-res-NetX2Y')
elseif Types==4 % coef perm
    Pvalu = ShowMatinfo.G2coefperm;
%     save test2
    Pval = Pvalu.P_mat1;
    BCCT_MatrixShowFinalResPermVer_other_GUI(Pval,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders,'GroupComp-CaSCN-Coef-X2Y')
    clear Pval
    Pval = Pvalu.P_mat2;
    BCCT_MatrixShowFinalResPermVer_other_GUI(Pval,Outputdir,enhanind,Coltype,PVALT,Pvalused,orders,'GroupComp-CaSCN-Coef-Y2X')
    clear Pval
end
end
function HSMATRIXEXIT(varargin)
Hshow = varargin{3};
close(Hshow.fig);
BCCT_VIEWmain;
end