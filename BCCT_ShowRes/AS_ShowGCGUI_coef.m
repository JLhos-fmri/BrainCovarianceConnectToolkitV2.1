function AS_ShowGCGUI_coef
[pat nam ext] = fileparts(which('AS_ShowFCGUI.m'));
templatepath = fullfile(pat,'mni_icbm152_t1_tal_nlin_asym_09a.nii');
Hsize = get(0,'screensize');
MIDPOINT = [Hsize(3)/2,Hsize(4)/2];
Asize = [100*3,100+40];
MaxSIZE = [Hsize(3) Hsize(4)]*0.8;
factor = MaxSIZE./Asize;
factornew = min(factor);
POSSIZE = Asize*factornew;
Hshow.fig = figure('position',[MIDPOINT(1)-POSSIZE(1)/2,MIDPOINT(2)-POSSIZE(2)/2,POSSIZE(1),POSSIZE(2)]);
Hshow.pat = pat;

Hshow.IO = uibuttongroup('parent',Hshow.fig,...
    'unit','norm',...
    'pos',[0 120/140 2/3 20/140]);
Hshow.Sagittal = uibuttongroup('parent',Hshow.fig,...
    'unit','norm',...
    'pos',[0 20/140 1/3 100/140]);
Hshow.Cornoral = uibuttongroup('parent',Hshow.fig,...
    'unit','norm',...
    'pos',[1/3 20/140 1/3 100/140]);
Hshow.Axial = uibuttongroup('parent',Hshow.fig,...
    'unit','norm',...
    'pos',[2/3 20/140 1/3 100/140]);
Hshow.Sagittal_Control = uibuttongroup('parent',Hshow.fig,...
    'unit','norm',...
    'pos',[0 0 1/3 20/140]);
Hshow.Cornoral_Control = uibuttongroup('parent',Hshow.fig,...
    'unit','norm',...
    'pos',[1/3 0 1/3 20/140]);
Hshow.Axial_Control = uibuttongroup('parent',Hshow.fig,...
    'unit','norm',...
    'pos',[2/3 0 1/3 20/140]);
Hshow.Other = uibuttongroup('parent',Hshow.fig,...
    'unit','norm',...
    'pos',[2/3 120/140 1/3 20/140]);
%%
% Hshow.IO
Hshow.IO_inputpb = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.05 0.725 0.1 0.2],'style','text','string','Overlay image');
Hshow.IO_inputsel = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.85 0.725 0.1 0.2],'style','pushbutton','string','Select');
Hshow.IO_inputed = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.2 0.725 0.6 0.2],'style','edit');

Hshow.IO_inputpb_p = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.05 0.5 0.1 0.2],'style','text','string','P image');
Hshow.IO_inputsel_p = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.85 0.5 0.1 0.2],'style','pushbutton','string','Select');
Hshow.IO_inputed_p = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.2 0.5 0.6 0.2],'style','edit','string','NULL');

Hshow.IO_tempputpb = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.05 0.25 0.1 0.2],'style','text','string','Template image');
Hshow.IO_tempputsel = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.85 0.25 0.1 0.2],'style','pushbutton','string','Select');
Hshow.IO_tempputed = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.2 0.25 0.6 0.2],'style','edit','string',templatepath);

Hshow.IO_Outputpb = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.05 0.025 0.1 0.2],'style','text','string','Outputdir');
Hshow.IO_Outputsel = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.85 0.025 0.1 0.2],'style','pushbutton','string','Select');
pwdpath = pwd;
Hshow.IO_Outputed = uicontrol('parent',Hshow.IO,'units','norm','pos',[0.2 0.025 0.6 0.2],'style','edit','string',pwdpath);

%%
Hshow.Other_show = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.075,0.1,0.2,0.3],'style','pushbutton','string','Show');
Hshow.Other_print = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.3,0.1,0.2,0.3],'style','pushbutton','string','Print All');
Hshow.Other_reset = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.525,0.1,0.2,0.3],'style','pushbutton','string','Reset');
Hshow.Other_exit = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.75,0.1,0.2,0.3],'style','pushbutton','string','Exit');

Hshow.Other_Threshold = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.075,0.8,0.2,0.15],'style','text','string','Threshold');
Hshow.Other_clus = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.3,0.8,0.2,0.15],'style','text','string','Cluster');
Hshow.Other_type = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.525,0.8,0.2,0.15],'style','text','string','GC');
Hshow.Other_dof = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.75,0.8,0.2,0.15],'style','text','string','Degree of Freedom');
set(Hshow.Other_type,'enable','off');
set(Hshow.Other_dof,'enable','off');
% Hshow.Other_Pth = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.075,0.4,0.1,0.15],'style','rad','string','P=','val',0);
Hshow.Other_Pth = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.075,0.4,0.1,0.15],'style','text','string','P=');
Hshow.Other_PthE = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.175,0.4,0.1,0.15],'style','edit');
Hshow.Other_PAlter = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.3,0.4,0.1,0.15],'style','pushbutton','string','Alter');
Hshow.Other_PFDR = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.425,0.4,0.1,0.15],'style','pushbutton','string','FDR');
set(Hshow.Other_Pth,'enable','off');
set(Hshow.Other_PthE,'enable','off');
% set(Hshow.Other_Pth,'val',0);
Hshow.Other_ShowTYPES = uibuttongroup('parent',Hshow.Other,'units','norm','pos',[0.55 0.4 0.4 0.15]);
Hshow.Other_ShowTYPE1 = uicontrol('parent',Hshow.Other_ShowTYPES,'units','norm','pos',[0.05,0.05,0.3,0.95],'style','rad','string','both','value',1);
Hshow.Other_ShowTYPE2 = uicontrol('parent',Hshow.Other_ShowTYPES,'units','norm','pos',[0.325,0.05,0.3,0.95],'style','rad','string','pos','value',0);
Hshow.Other_ShowTYPE3 = uicontrol('parent',Hshow.Other_ShowTYPES,'units','norm','pos',[0.65,0.05,0.3,0.95],'style','rad','string','neg','value',0);

Hshow.Other_Threshold1 = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.075,0.6,0.1,0.2],'style','edit');
Hshow.Other_Threshold2 = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.175,0.6,0.1,0.2],'style','edit');
Hshow.Other_Cnum = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.3,0.6,0.2,0.2],'style','edit','string','5');
Hshow.Other_typeE = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.525,0.6,0.2,0.2],'style','edit','string','GC');
Hshow.Other_dofE = uicontrol('parent',Hshow.Other,'units','norm','pos',[0.75,0.6,0.2,0.2],'style','edit','string','Degree of Freedom','string','10000');
set(Hshow.Other_typeE,'enable','off');
set(Hshow.Other_dofE,'enable','off');
%%
Hshow.Sagittal_axes = axes('parent',Hshow.Sagittal,'units','norm','pos',[0.05 0.05 0.9 0.9]);
Hshow.Cornoral_axes = axes('parent',Hshow.Cornoral,'units','norm','pos',[0.05 0.05 0.9 0.9]);
Hshow.Axial_axes = axes('parent',Hshow.Axial,'units','norm','pos',[0.05 0.05 0.9 0.9]);

%%
Hshow.Sagittal_contslide = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.1 0.85 0.8 0.15],'style','slider','Min',1,'Max',197,'value',99,'SliderStep',[1/196 1]);
Hshow.Cornoral_contslide = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.1 0.85 0.8 0.15],'style','slider','Min',1,'Max',233,'value',135,'SliderStep',[1/232 1]);
Hshow.Axial_contslide = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.1 0.85 0.8 0.15],'style','slider','Min',1,'Max',189,'value',73,'SliderStep',[1/188  1]);
Hshow.Sagittal_mnitext = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.05 0.6 0.2 0.2],'style','text','string','mni: x=0');
Hshow.Cornoral_mnitext = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.05 0.6 0.2 0.2],'style','text','string','mni: y=0');
Hshow.Axial_mnitext = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.05 0.6 0.2 0.2],'style','text','string','mni: z=0');
Hshow.Sagittal_mniedit = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.25 0.6 0.1 0.2],'style','edit','string','0');
Hshow.Cornoral_mniedit = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.25 0.6 0.1 0.2],'style','edit','string','0');
Hshow.Axial_mniedit = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.25 0.6 0.1 0.2],'style','edit','string','0');

Hshow.Sagittal_cortext = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.4 0.6 0.2 0.2],'style','text','string','cor: x=99');
Hshow.Cornoral_cortext = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.4 0.6 0.2 0.2],'style','text','string','cor: y=135');
Hshow.Axial_cortext = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.4 0.6 0.2 0.2],'style','text','string','cor: z=73');
Hshow.Sagittal_coredit = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.65 0.6 0.1 0.2],'style','edit','string','99');
Hshow.Cornoral_coredit = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.65 0.6 0.1 0.2],'style','edit','string','135');
Hshow.Axial_coredit = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.65 0.6 0.1 0.2],'style','edit','string','73');

Hshow.Sagittal_change = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.85 0.6 0.1 0.2],'style','pushbutton','string','change');
Hshow.Cornoral_change = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.85 0.6 0.1 0.2],'style','pushbutton','string','change');
Hshow.Axial_change = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.85 0.6 0.1 0.2],'style','pushbutton','string','change');


Hshow.Sagittal_pa = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.05 0.2 0.2 0.2],'style','pushbutton','string','print all');
Hshow.Cornoral_pa = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.05 0.2 0.2 0.2],'style','pushbutton','string','print all');
Hshow.Axial_pa = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.05 0.2 0.2 0.2],'style','pushbutton','string','print all');
Hshow.Sagittal_pc = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.25 0.2 0.2 0.2],'style','pushbutton','string','print current');
Hshow.Cornoral_pc = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.25 0.2 0.2 0.2],'style','pushbutton','string','print current');
Hshow.Axial_pc = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.25 0.2 0.2 0.2],'style','pushbutton','string','print current');


Hshow.Sagittal_setrag = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.55 0.3 0.2 0.15],'style','pushbutton','string','Selrange:Cor');
Hshow.Cornoral_setrag = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.55 0.3 0.2 0.15],'style','pushbutton','string','Selrange:Cor');
Hshow.Axial_setrag = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.55 0.3 0.2 0.15],'style','pushbutton','string','Selrange:Cor');

Hshow.Sagittal_edrag = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.55 0.1 0.2 0.2],'style','edit');
Hshow.Cornoral_edrag = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.55 0.1 0.2 0.2],'style','edit');
Hshow.Axial_edrag = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.55 0.1 0.2 0.2],'style','edit');

Hshow.Sagittal_prrag = uicontrol('parent',Hshow.Sagittal_Control,'units','norm','pos',[0.75 0.1 0.2 0.35],'style','pushbutton','string','print');
Hshow.Cornoral_prrag = uicontrol('parent',Hshow.Cornoral_Control,'units','norm','pos',[0.75 0.1 0.2 0.35],'style','pushbutton','string','print');
Hshow.Axial_prrag = uicontrol('parent',Hshow.Axial_Control,'units','norm','pos',[0.75 0.1 0.2 0.35],'style','pushbutton','string','print');

load(fullfile(pat,'DefaultTemplateMatForShow','Vbg.mat'));
Tempshowbg = load(fullfile(pat,'DefaultTemplateMatForShow','Z73.mat'));
AxialDAT = rot90(Tempshowbg.Slice);
sizeaxial = size(AxialDAT);
difaxis = sizeaxial(1)-sizeaxial(2);
AxialDATnew = ones(max(sizeaxial))*30;
if difaxis>0
    midc = floor(difaxis/2);
    AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
else
    midc = floor(-difaxis/2);
    AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
end
%
Tempshowbg = load(fullfile(pat,'DefaultTemplateMatForShow','Y135.mat'));
% AxialDAT = rot90(Tempshowbg.Slice);
% CornoralDAT = rot90(squeeze(DbgRe(:,135,:)));
CornoralDAT = rot90(Tempshowbg.Slice);
sizeCornoral = size(CornoralDAT);
difcornoral = sizeCornoral(1)-sizeCornoral(2);
CornoralDATnew = ones(max(sizeCornoral))*30;
if difcornoral>0
    midc = floor(difcornoral/2);
    CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
else
    midc = floor(-difcornoral/2);
    CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
end
%
Tempshowbg = load(fullfile(pat,'DefaultTemplateMatForShow','X99.mat'));
SagittalDAT = rot90(Tempshowbg.Slice);
% SagittalDAT = rot90(squeeze(DbgRe(99,:,:)));
sizeSagittal = size(SagittalDAT);
difsagittal = sizeSagittal(1)-sizeSagittal(2);
SagittalDATnew = ones(max(sizeSagittal))*30;
if difsagittal>0
    midc = floor(difsagittal/2);
    SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
else
    midc = floor(-difsagittal/2);
    SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
end
% Hshow.DbgRe = DbgRe;
Hshow.Vbg = Vbg;

if isempty(dir([Hshow.pat,filesep,'TempForOrigShow']))
    mkdir([Hshow.pat,filesep,'TempForOrigShow']);
else
    rmdir([Hshow.pat,filesep,'TempForOrigShow'],'s');
    mkdir([Hshow.pat,filesep,'TempForOrigShow']);
end

copyfile([pat,filesep,'DefaultTemplateMatForShow',filesep,'*.mat'],[pat,filesep,'TempForOrigShow']);

% image(backbin2,'parent',ASBC.mainaxes,'CDataMapping','scaled');
image(AxialDATnew,'parent',Hshow.Axial_axes,'CDataMapping','scaled');axis(Hshow.Axial_axes,'off');colormap(Hshow.Axial_axes,'gray')
image(CornoralDATnew,'parent',Hshow.Cornoral_axes,'CDataMapping','scaled');axis(Hshow.Cornoral_axes,'off');colormap(Hshow.Cornoral_axes,'gray')
image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');colormap(Hshow.Sagittal_axes,'gray')
text(size(AxialDATnew,1),size(AxialDATnew,1)/2,'R','parent',Hshow.Axial_axes)
text(size(CornoralDATnew,1),size(CornoralDATnew,1)/2,'R','parent',Hshow.Cornoral_axes)
text(size(SagittalDATnew,1),size(SagittalDATnew,1)/2,'A','parent',Hshow.Sagittal_axes)
set(Hshow.Axial_axes,'Clim',[30 85]);
set(Hshow.Cornoral_axes,'Clim',[30 85]);
set(Hshow.Sagittal_axes,'Clim',[30 85]);

%% add the callback
set(Hshow.IO_inputsel,'callback',{@ASmapshowInput,Hshow});
set(Hshow.IO_inputsel_p,'callback',{@ASmapshowInputP,Hshow});
set(Hshow.IO_tempputsel,'callback',{@ASmapshowTemplate,Hshow});
set(Hshow.IO_Outputsel,'callback',{@ASmapshowoutput,Hshow});
set(Hshow.IO_inputed,'callback',{@ASmapshowinputedit,Hshow});
set(Hshow.IO_inputed_p,'callback',{@ASmapshowinputeditp,Hshow});
set(Hshow.IO_tempputed,'callback',{@ASmapshowTemplateedit,Hshow});
%
set(Hshow.Other_PFDR,'callback',{@ASmapshowfdr,Hshow});
set(Hshow.Other_show,'callback',{@ASmapshowshow,Hshow});
set(Hshow.Other_print,'callback',{@ASmapshowprint,Hshow});
set(Hshow.Other_reset,'callback',{@ASmapshowreset,Hshow});
set(Hshow.Other_exit,'callback',{@ASmapshowexit,Hshow});
set(Hshow.Other_PAlter,'callback',{@ASmapshowPvalchange,Hshow});
%
set(Hshow.Sagittal_contslide,'callback',{@ASmapshowSagConSlid,Hshow});
set(Hshow.Cornoral_contslide,'callback',{@ASmapshowCorConSlid,Hshow});
set(Hshow.Axial_contslide,'callback',{@ASmapshowAxiConSlid,Hshow});
set(Hshow.Sagittal_change,'callback',{@ASmapshowSagChange,Hshow});
set(Hshow.Cornoral_change,'callback',{@ASmapshowCorChange,Hshow});
set(Hshow.Axial_change,'callback',{@ASmapshowAxiChange,Hshow});

set(Hshow.Sagittal_pa,'callback',{@ASmapshowSagPrintAll,Hshow});
set(Hshow.Cornoral_pa,'callback',{@ASmapshowCorPrintAll,Hshow});
set(Hshow.Axial_pa,'callback',{@ASmapshowAxiPrintAll,Hshow});
set(Hshow.Sagittal_pc,'callback',{@ASmapshowSagPrintCurrent,Hshow});
set(Hshow.Cornoral_pc,'callback',{@ASmapshowCorPrintCurrent,Hshow});
set(Hshow.Axial_pc,'callback',{@ASmapshowAxiPrintCurrent,Hshow});

set(Hshow.Sagittal_setrag,'callback',{@ASmapshowSagRangeshow,Hshow});
set(Hshow.Cornoral_setrag,'callback',{@ASmapshowCorRangeshow,Hshow});
set(Hshow.Axial_setrag,'callback',{@ASmapshowAxiRangeshow,Hshow});
set(Hshow.Sagittal_edrag,'callback',{@ASmapshowSagRag_ed,Hshow});
set(Hshow.Cornoral_edrag,'callback',{@ASmapshowCorRag_ed,Hshow});
set(Hshow.Axial_edrag,'callback',{@ASmapshowAxiRag_ed,Hshow});

set(Hshow.Sagittal_prrag,'callback',{@ASmapshowSagPrintRange,Hshow});
set(Hshow.Cornoral_prrag,'callback',{@ASmapshowCorPrintRange,Hshow});
set(Hshow.Axial_prrag,'callback',{@ASmapshowAxiPrintRange,Hshow});
end
%%
function ASmapshowInput(varargin)
Hshow = varargin{3};
[nam,path,ext] = uigetfile({'*.nii';'*.img';'*.hdr';'*.*'},'Image directory for show');
patnam = fullfile(path,nam);
set(Hshow.IO_inputed,'string',patnam);
% AS_show_gainIO(Hshow);
pat = Hshow.pat;
% load(fullfile(pat,'InfoofStatForShow','DOFinfo.mat'))
% set(Hshow.Other_typeE,'string',styleofstat);
% DOF = [dof1,dof2];
% set(Hshow.Other_dofE,'string',num2str(DOF));
Pthv = 0.05;
set(Hshow.Other_PthE,'string','0.05');
% [Zout,out] = AS_PtoTRF(Pthv,styleofstat,dof1,dof2);
indir = get(Hshow.IO_inputed,'string');
[vmap,dmap] = Dynamic_read_dir_NIFTI(indir);

dmap(isnan(dmap)) = 0;
dmap(isinf(dmap)) = 0;
maxv = max(abs(dmap));
out = maxv/2;
save('test.mat')
set(Hshow.Other_Threshold1,'string',num2str(out));
set(Hshow.Other_Threshold2,'string',num2str(maxv));

% if strcmpi(styleofstat,'P')
%     set(Hshow.Other_Threshold1,'string','0.95');
%     set(Hshow.Other_Threshold2,'string','1');
% else    
%     dmap(isnan(dmap)) = 0;
%     dmap(isinf(dmap)) = 0;
%     maxv = max(abs(dmap));
%     set(Hshow.Other_Threshold1,'string',num2str(out));
%     set(Hshow.Other_Threshold2,'string',num2str(maxv));
% end
end
function ASmapshowinputedit(varargin)
Hshow = varargin{3};
% AS_show_gainIO(Hshow);
% pat = Hshow.pat;
% load(fullfile(pat,'InfoofStatForShow','DOFinfo.mat'))
% set(Hshow.Other_typeE,'string',styleofstat);
% DOF = [dof1,dof2];
% set(Hshow.Other_dofE,'string',num2str(DOF));
% Pthv = 0.05;
set(Hshow.Other_PthE,'string','0.05');

indir = get(Hshow.IO_inputed,'string');
[vmap,dmap] = Dynamic_read_dir_NIFTI(indir);
dmap(isnan(dmap)) = 0;
dmap(isinf(dmap)) = 0;
maxv = max(abs(dmap));
out = maxv/2;
set(Hshow.Other_Threshold1,'string',num2str(out));
set(Hshow.Other_Threshold2,'string',num2str(maxv));

end
function ASmapshowInputP(varargin)
Hshow = varargin{3};
[nam,path,ext] = uigetfile({'*.nii';'*.img';'*.hdr';'*.*'},'Image directory for show');
patnam = fullfile(path,nam);
set(Hshow.IO_inputed_p,'string',patnam);

set(Hshow.Other_Pth,'enable','on');
set(Hshow.Other_PthE,'enable','on');

indir = get(Hshow.IO_inputed,'string');
[vmap,dmap] = Dynamic_read_dir_NIFTI(indir);
[vmap_p,dmap_p] = Dynamic_read_dir_NIFTI(patnam);
pthr = 0.05;
save('testp.mat')
out = min(dmap(dmap_p>=(1-pthr)&dmap_p>0&abs(dmap)>0));
maxv = max(abs(dmap));
set(Hshow.Other_Threshold1,'string',num2str(out));
set(Hshow.Other_Threshold2,'string',num2str(maxv));

end
function ASmapshowInputeditp(varargin)
Hshow = varargin{3};
set(Hshow.Other_Pth,'enable','on');
set(Hshow.Other_PthE,'enable','on');
patnam = get(Hshow.IO_inputed_p,'string');
indir = get(Hshow.IO_inputed,'string');
[vmap,dmap] = Dynamic_read_dir_NIFTI(indir);
[vmap_p,dmap_p] = Dynamic_read_dir_NIFTI(patnam);
pthr = 0.05;
dmapu = dmap(dmap_p<=pthr|dmap_p>=1-pthr);
out = min(abs(dmapu));
% out = min(dmap(dmap_p>=(1-pthr)));
maxv = max(abs(dmap));
set(Hshow.Other_Threshold1,'string',num2str(out));
set(Hshow.Other_Threshold2,'string',num2str(maxv));
end

function ASmapshowTemplateedit(varargin)
Hshow = varargin{3};
% [nam,path,ext] = uigetfile({'*.nii';'*.img';'*.hdr';'*.*'},'template directory for show');
% patnam = get(Hshow.IO_tempputed,'string');
% set(Hshow.IO_tempputed,'string',patnam);
AS_show_gainIO_temp(Hshow);

load(fullfile(pat,'TempForOrigShow','Vbg.mat'));
mnicoord = [0 0 0];
corcoord = mni2cor(mnicoord,Vbg.mat);

Tempshowbg = load(fullfile(pat,'TempForOrigShow',['Z',num2str(corcoord(3)),'.mat']));
AxialDAT = rot90(Tempshowbg.Slice);
sizeaxial = size(AxialDAT);
difaxis = sizeaxial(1)-sizeaxial(2);
AxialDATnew = ones(max(sizeaxial))*30;
if difaxis>0
    midc = floor(difaxis/2);
    AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
else
    midc = floor(-difaxis/2);
    AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
end
%
Tempshowbg = load(fullfile(pat,'TempForOrigShow','Y',[num2str(corcoord(32)),'.mat']));
% AxialDAT = rot90(Tempshowbg.Slice);
% CornoralDAT = rot90(squeeze(DbgRe(:,135,:)));
CornoralDAT = rot90(Tempshowbg.Slice);
sizeCornoral = size(CornoralDAT);
difcornoral = sizeCornoral(1)-sizeCornoral(2);
CornoralDATnew = ones(max(sizeCornoral))*30;
if difcornoral>0
    midc = floor(difcornoral/2);
    CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
else
    midc = floor(-difcornoral/2);
    CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
end
%
Tempshowbg = load(fullfile(pat,'TempForOrigShow',['X',num2str(corcoord(1)),'.mat']));
SagittalDAT = rot90(Tempshowbg.Slice);
% SagittalDAT = rot90(squeeze(DbgRe(99,:,:)));
sizeSagittal = size(SagittalDAT);
difsagittal = sizeSagittal(1)-sizeSagittal(2);
SagittalDATnew = ones(max(sizeSagittal))*30;
if difsagittal>0
    midc = floor(difsagittal/2);
    SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
else
    midc = floor(-difsagittal/2);
    SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
end
end
function ASmapshowTemplate(varargin)
Hshow = varargin{3};
[nam,path,ext] = uigetfile({'*.nii';'*.img';'*.hdr';'*.*'},'template directory for show');
patnam = fullfile(path,nam);
set(Hshow.IO_tempputed,'string',patnam);
AS_show_gainIO_temp(Hshow);

load(fullfile(pat,'TempForOrigShow','Vbg.mat'));
mnicoord = [0 0 0];
corcoord = mni2cor(mnicoord,Vbg.mat);

Tempshowbg = load(fullfile(pat,'TempForOrigShow',['Z',num2str(corcoord(3)),'.mat']));
AxialDAT = rot90(Tempshowbg.Slice);
sizeaxial = size(AxialDAT);
difaxis = sizeaxial(1)-sizeaxial(2);
AxialDATnew = ones(max(sizeaxial))*30;
if difaxis>0
    midc = floor(difaxis/2);
    AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
else
    midc = floor(-difaxis/2);
    AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
end
%
Tempshowbg = load(fullfile(pat,'TempForOrigShow','Y',[num2str(corcoord(32)),'.mat']));
% AxialDAT = rot90(Tempshowbg.Slice);
% CornoralDAT = rot90(squeeze(DbgRe(:,135,:)));
CornoralDAT = rot90(Tempshowbg.Slice);
sizeCornoral = size(CornoralDAT);
difcornoral = sizeCornoral(1)-sizeCornoral(2);
CornoralDATnew = ones(max(sizeCornoral))*30;
if difcornoral>0
    midc = floor(difcornoral/2);
    CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
else
    midc = floor(-difcornoral/2);
    CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
end
%
Tempshowbg = load(fullfile(pat,'TempForOrigShow',['X',num2str(corcoord(1)),'.mat']));
SagittalDAT = rot90(Tempshowbg.Slice);
% SagittalDAT = rot90(squeeze(DbgRe(99,:,:)));
sizeSagittal = size(SagittalDAT);
difsagittal = sizeSagittal(1)-sizeSagittal(2);
SagittalDATnew = ones(max(sizeSagittal))*30;
if difsagittal>0
    midc = floor(difsagittal/2);
    SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
else
    midc = floor(-difsagittal/2);
    SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
end
end
function ASmapshowoutput(varargin)
Hshow = varargin{3};
Path = uigetdir(pwd,'Output for pic and other information');
set(Hshow.IO_Outputed,'string',Path);
end
%%
function ASmapshowfdr(varargin)
Hshow = varargin{3};
uiwait(msgbox('comming soon'));
end
function ASmapshowshow(varargin)
Hshow = varargin{3};
inputdir = get(Hshow.IO_inputed,'string');
inputdirP = get(Hshow.IO_inputed_p,'string');
templatedir = get(Hshow.IO_tempputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');
Pval = get(Hshow.Other_PthE,'string');
Types = get(Hshow.Other_typeE,'string');
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Clusnum = get(Hshow.Other_Cnum,'string');
SliceSagX = get(Hshow.Sagittal_contslide,'val');
SliceCorX = get(Hshow.Cornoral_contslide,'val');
SliceAxiX = get(Hshow.Axial_contslide,'val');

SliceSag = uint8(SliceSagX);
SliceCor = uint8(SliceCorX);
SliceAxi = uint8(SliceAxiX);

% DOF = get(Hshow.Other_dofE,'string');

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);

pthrv = str2num(Pval);
colmaxt = [str2num(Thrval1),str2num(Thrval2)];
CLUSNUM = str2num(Clusnum);
% [DOUTSHOWaxial,DOUTSHOWcornoral,DOUTSHOWsagittal,dshownum] = AS_MapFCshow3sliceUP_GC_res(outdir,inputdir,inputdirP,templatedir,pthrv,colmaxt,CLUSNUM);
[DOUTSHOWaxial,DOUTSHOWcornoral,DOUTSHOWsagittal,DOUTSHOWaxialPOS,DOUTSHOWcornoralPOS,...
    DOUTSHOWsagittalPOS,DOUTSHOWaxialNEG,DOUTSHOWcornoralNEG,DOUTSHOWsagittalNEG,...
    dshownum] = AS_MapFCshow3sliceUP_GC_coef(outdir,inputdir,inputdirP,templatedir,pthrv,colmaxt,CLUSNUM);

if Showtype1==1
    SagShowmap = DOUTSHOWsagittal(:,:,SliceSag);
    CorShowmap = DOUTSHOWcornoral(:,:,SliceCor);
    AxiShowmap = DOUTSHOWaxial(:,:,SliceAxi);
    clear DOUTSHOWaxial DOUTSHOWcornoral DOUTSHOWsagittal DOUTSHOWaxialPOS DOUTSHOWcornoralPOS DOUTSHOWsagittalPOS DOUTSHOWaxialNEG DOUTSHOWcornoralNEG DOUTSHOWsagittalNEG
    AxialDAT = AxiShowmap;
    sizeaxial = size(AxialDAT);
    difaxis = sizeaxial(1)-sizeaxial(2);
    AxialDATnew = ones(max(sizeaxial))*(-dshownum);
    if difaxis>0
        midc = floor(difaxis/2);
        AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
    else
        midc = floor(-difaxis/2);
        AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
    end
    %
    CornoralDAT = CorShowmap;
    sizeCornoral = size(CornoralDAT);
    difcornoral = sizeCornoral(1)-sizeCornoral(2);
    CornoralDATnew = ones(max(sizeCornoral))*(-dshownum);
    if difcornoral>0
        midc = floor(difcornoral/2);
        CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
    else
        midc = floor(-difcornoral/2);
        CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
    end
    %
    SagittalDAT = SagShowmap;
    sizeSagittal = size(SagittalDAT);
    difsagittal = sizeSagittal(1)-sizeSagittal(2);
    SagittalDATnew = ones(max(sizeSagittal))*(-dshownum);
    if difsagittal>0
        midc = floor(difsagittal/2);
        SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
    else
        midc = floor(-difsagittal/2);
        SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
    end
    % image(backbin2,'parent',ASBC.mainaxes,'CDataMapping','scaled');
    image(AxialDATnew,'parent',Hshow.Axial_axes,'CDataMapping','scaled');axis(Hshow.Axial_axes,'off');colormap(Hshow.Axial_axes,COLMAPBOTH)
    image(CornoralDATnew,'parent',Hshow.Cornoral_axes,'CDataMapping','scaled');axis(Hshow.Cornoral_axes,'off');colormap(Hshow.Axial_axes,COLMAPBOTH)
    image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');colormap(Hshow.Axial_axes,COLMAPBOTH)
    set(Hshow.Axial_axes,'Clim',[-dshownum*2 dshownum*2]);
    set(Hshow.Cornoral_axes,'Clim',[-dshownum*2 dshownum*2]);
    set(Hshow.Sagittal_axes,'Clim',[-dshownum*2 dshownum*2]);
elseif Showtype2==1 %pos only
    
    SagShowmap = DOUTSHOWsagittalPOS(:,:,SliceSag);
    CorShowmap = DOUTSHOWcornoralPOS(:,:,SliceCor);
    AxiShowmap = DOUTSHOWaxialPOS(:,:,SliceAxi);
    clear DOUTSHOWaxial DOUTSHOWcornoral DOUTSHOWsagittal DOUTSHOWaxialPOS DOUTSHOWcornoralPOS DOUTSHOWsagittalPOS DOUTSHOWaxialNEG DOUTSHOWcornoralNEG DOUTSHOWsagittalNEG
    AxialDAT = AxiShowmap;
    %         AxialDAT(AxiShowmap<-dshownum) = -dshownum;
    sizeaxial = size(AxialDAT);
    difaxis = sizeaxial(1)-sizeaxial(2);
    AxialDATnew = ones(max(sizeaxial))*(-dshownum);
    if difaxis>0
        midc = floor(difaxis/2);
        AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
    else
        midc = floor(-difaxis/2);
        AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
    end
    %
    CornoralDAT = CorShowmap;
    %         CornoralDAT(CorShowmap<-dshownum) = -dshownum;
    sizeCornoral = size(CornoralDAT);
    difcornoral = sizeCornoral(1)-sizeCornoral(2);
    CornoralDATnew = ones(max(sizeCornoral))*(-dshownum);
    if difcornoral>0
        midc = floor(difcornoral/2);
        CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
    else
        midc = floor(-difcornoral/2);
        CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
    end
    %
    SagittalDAT = SagShowmap;
    %         SagittalDAT(SagShowmap<-dshownum) = -dshownum;
    sizeSagittal = size(SagittalDAT);
    difsagittal = sizeSagittal(1)-sizeSagittal(2);
    SagittalDATnew = ones(max(sizeSagittal))*(-dshownum);
    if difsagittal>0
        midc = floor(difsagittal/2);
        SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
    else
        midc = floor(-difsagittal/2);
        SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
    end
    % image(backbin2,'parent',ASBC.mainaxes,'CDataMapping','scaled');
    image(AxialDATnew,'parent',Hshow.Axial_axes,'CDataMapping','scaled');axis(Hshow.Axial_axes,'off');colormap(Hshow.Axial_axes,COLMAPPOS)
    image(CornoralDATnew,'parent',Hshow.Cornoral_axes,'CDataMapping','scaled');axis(Hshow.Cornoral_axes,'off');colormap(Hshow.Axial_axes,COLMAPPOS)
    image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');colormap(Hshow.Axial_axes,COLMAPPOS)
    set(Hshow.Axial_axes,'Clim',[-dshownum*1 dshownum*2]);
    set(Hshow.Cornoral_axes,'Clim',[-dshownum*1 dshownum*2]);
    set(Hshow.Sagittal_axes,'Clim',[-dshownum*1 dshownum*2]);
    
elseif Showtype3==1
    SagShowmap = DOUTSHOWsagittalNEG(:,:,SliceSag);
    CorShowmap = DOUTSHOWcornoralNEG(:,:,SliceCor);
    AxiShowmap = DOUTSHOWaxialNEG(:,:,SliceAxi);
    clear DOUTSHOWaxial DOUTSHOWcornoral DOUTSHOWsagittal DOUTSHOWaxialPOS DOUTSHOWcornoralPOS DOUTSHOWsagittalPOS DOUTSHOWaxialNEG DOUTSHOWcornoralNEG DOUTSHOWsagittalNEG
    AxialDAT = AxiShowmap;
    %         AxialDAT(AxiShowmap>dshownum) = dshownum;
    sizeaxial = size(AxialDAT);
    difaxis = sizeaxial(1)-sizeaxial(2);
    AxialDATnew = ones(max(sizeaxial))*(-dshownum);
    if difaxis>0
        midc = floor(difaxis/2);
        AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
    else
        midc = floor(-difaxis/2);
        AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
    end
    %
    CornoralDAT = CorShowmap;
    %         CornoralDAT(CorShowmap>dshownum) = dshownum;
    sizeCornoral = size(CornoralDAT);
    difcornoral = sizeCornoral(1)-sizeCornoral(2);
    CornoralDATnew = ones(max(sizeCornoral))*(-dshownum);
    if difcornoral>0
        midc = floor(difcornoral/2);
        CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
    else
        midc = floor(-difcornoral/2);
        CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
    end
    %
    SagittalDAT = SagShowmap;
    %         SagittalDAT(SagShowmap>dshownum) = dshownum;
    sizeSagittal = size(SagittalDAT);
    difsagittal = sizeSagittal(1)-sizeSagittal(2);
    SagittalDATnew = ones(max(sizeSagittal))*(-dshownum);
    if difsagittal>0
        midc = floor(difsagittal/2);
        SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
    else
        midc = floor(-difsagittal/2);
        SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
    end
    % image(backbin2,'parent',ASBC.mainaxes,'CDataMapping','scaled');
    image(AxialDATnew,'parent',Hshow.Axial_axes,'CDataMapping','scaled');axis(Hshow.Axial_axes,'off');colormap(Hshow.Axial_axes,COLMAPNEG)
    image(CornoralDATnew,'parent',Hshow.Cornoral_axes,'CDataMapping','scaled');axis(Hshow.Cornoral_axes,'off');colormap(Hshow.Axial_axes,COLMAPNEG)
    image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');colormap(Hshow.Axial_axes,COLMAPNEG)
    set(Hshow.Axial_axes,'Clim',[-dshownum*2 dshownum*1]);
    set(Hshow.Cornoral_axes,'Clim',[-dshownum*2 dshownum*1]);
    set(Hshow.Sagittal_axes,'Clim',[-dshownum*2 dshownum*1]);
end

text(size(AxialDATnew,1),size(AxialDATnew,1)/2,'R','parent',Hshow.Axial_axes)
text(size(CornoralDATnew,1),size(CornoralDATnew,1)/2,'R','parent',Hshow.Cornoral_axes)
text(size(SagittalDATnew,1),size(SagittalDATnew,1)/2,'A','parent',Hshow.Sagittal_axes)
uiwait(msgbox('load OK'));
end
function ASmapshowprint(varargin)
Hshow = varargin{3};
uiwait(msgbox('We will print every slice of sagittal/Cornoral/Axial images'));

inputdir = get(Hshow.IO_inputed,'string');
templatedir = get(Hshow.IO_tempputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');
Pval = get(Hshow.Other_PthE,'string');
Types = get(Hshow.Other_typeE,'string');
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Clusnum = get(Hshow.Other_Cnum,'string');
SliceSag = get(Hshow.Sagittal_contslide,'val');
SliceCor = get(Hshow.Cornoral_contslide,'val');
SliceAxi = get(Hshow.Axial_contslide,'val');
DOF = get(Hshow.Other_dofE,'string');

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);

DOFv = str2num(DOF);
if strcmpi(Types,'F')&&length(DOFv)<2
    error('wrong Degree of Freedom for F stastical analysis');
end
dof1 = DOFv(1);
if length(DOFv)<2
    dof20 = [];
else
    dof20 = DOFv(2);
end
pthrv = str2num(Pval);
colmaxt = [str2num(Thrval1),str2num(Thrval2)];
CLUSNUM = str2num(Clusnum);
% [DOUTSHOWaxial,DOUTSHOWcornoral,DOUTSHOWsagittal,dshownum]=AS_MapFCshow3sliceU(outdir,inputdir,templatedir,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM);
% SagShowmap = DOUTSHOWsagittal(:,:,SliceSag);
% CorShowmap = DOUTSHOWcornoral(:,:,SliceCor);
% AxiShowmap = DOUTSHOWaxial(:,:,SliceAxi);
if strcmpi(Types,'P')
    ShowtypeTry = [Showtype1,Showtype2,Showtype3];
    OutnameLab = 'PrintAllSlice';
    AS_MapFCshow_used_simpleUP(outdir,inputdir,templatedir,OutnameLab,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM,ShowtypeTry)
    
else
    ShowtypeTry = [Showtype1,Showtype2,Showtype3];
    OutnameLab = 'PrintAllSlice';
    AS_MapFCshow_used_simpleU(outdir,inputdir,templatedir,OutnameLab,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM,ShowtypeTry)
    % AS_MapFCshow_used_simple(outdir,inputdir,OutnameLab,pthrv,labmarkinput,dof1,dof20,colmaxt,CLUSNUM)
end
end
function ASmapshowreset(varargin)
Hshow = varargin{3};
close(Hshow.fig);
AS_ShowFCGUI
end
function ASmapshowexit(varargin)
Hshow = varargin{3};
close(Hshow.fig);
ASBC_VIEW;
end
%%
function ASmapshowPvalchange(varargin)
Hshow = varargin{3};
% DOF = get(Hshow.Other_dofE,'string');
% DOFv = str2num(DOF);
% Types = get(Hshow.Other_typeE,'string');
save ttt1
Pval = get(Hshow.Other_PthE,'string');
pthr = str2num(Pval);
patnam = get(Hshow.IO_inputed_p,'string');
indir = get(Hshow.IO_inputed,'string');
[vmap,dmap] = Dynamic_read_dir_NIFTI(indir);
[vmap_p,dmap_p] = Dynamic_read_dir_NIFTI(patnam);
% pthr = 0.05;
out = min(dmap(dmap_p>=(1-pthr)));
maxv = max(abs(dmap));
set(Hshow.Other_Threshold1,'string',num2str(out));
set(Hshow.Other_Threshold2,'string',num2str(maxv));
end
function ASmapshowSagConSlid(varargin)
Hshow = varargin{3};
inputdir = get(Hshow.IO_inputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
SagVal = get(Hshow.Sagittal_contslide,'val');
SagVal = uint8(SagVal);
corpos = [double(SagVal),1.0,1.0];
mnipos = cor2mni(corpos,Hshow.Vbg.mat);
MNI = mnipos(1);
set(Hshow.Sagittal_coredit,'string',num2str(SagVal));
set(Hshow.Sagittal_mniedit,'string',num2str(MNI));
set(Hshow.Sagittal_cortext,'string',['cor:x=',num2str(SagVal)]);
set(Hshow.Sagittal_mnitext,'string',['mni:x=',num2str(MNI)]);

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');


if isempty(inputdir)
    load([Hshow.pat,filesep,'TempForOrigShow',filesep,'X',num2str(SagVal),'.mat']);
    SagittalDAT = rot90(Slice);
    sizeSagittal = size(SagittalDAT);
    difsagittal = sizeSagittal(1)-sizeSagittal(2);
    SagittalDATnew = ones(max(sizeSagittal))*30;
    if difsagittal>0
        midc = floor(difsagittal/2);
        SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
    else
        midc = floor(-difsagittal/2);
        SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
    end
    image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');colormap(Hshow.Sagittal_axes,'gray')
    set(Hshow.Sagittal_axes,'Clim',[30 85]);
else
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(1)-colmaxt(2));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        SagittalDAT = DOUTSHOWsagittal(:,:,SagVal);
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        SagittalDAT = DOUTSHOWsagittalPOS(:,:,SagVal);
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        SagittalDAT = DOUTSHOWsagittalNEG(:,:,SagVal);
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG
    end
    sizeSagittal = size(SagittalDAT);
    difsagittal = sizeSagittal(1)-sizeSagittal(2);
    SagittalDATnew = ones(max(sizeSagittal))*-dshownum;
    if difsagittal>0
        midc = floor(difsagittal/2);
        SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
    else
        midc = floor(-difsagittal/2);
        SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
    end
    image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');
    if Showtype1
        set(Hshow.Sagittal_axes,'Clim',[-dshownum*2 dshownum*2]);
        colormap(COLMAPBOTH)
    elseif Showtype2
        set(Hshow.Sagittal_axes,'Clim',[-dshownum*1 dshownum*2]);
        colormap(COLMAPPOS)        
    elseif Showtype3
        set(Hshow.Sagittal_axes,'Clim',[-dshownum*2 dshownum*1]);
        colormap(COLMAPNEG)        
    end    
end
end
function ASmapshowCorConSlid(varargin)
Hshow = varargin{3};
inputdir = get(Hshow.IO_inputed,'string');
Corval = get(Hshow.Cornoral_contslide,'val');
outdir = get(Hshow.IO_Outputed,'string');
Corval = uint8(Corval);
corpos = [1.0,double(Corval),1.0];
mnipos = cor2mni(corpos,Hshow.Vbg.mat);
MNI = mnipos(2);
set(Hshow.Cornoral_coredit,'string',num2str(Corval));
set(Hshow.Cornoral_mniedit,'string',num2str(MNI));
set(Hshow.Cornoral_cortext,'string',['cor:y=',num2str(Corval)]);
set(Hshow.Cornoral_mnitext,'string',['mni:y=',num2str(MNI)]);

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');


if isempty(inputdir)
    load([Hshow.pat,filesep,'TempForOrigShow',filesep,'Y',num2str(Corval),'.mat']);
    CornoralDAT = rot90(Slice);
    sizeCornoral = size(CornoralDAT);
    difcornoral = sizeCornoral(1)-sizeCornoral(2);
    CornoralDATnew = ones(max(sizeCornoral))*30;
    if difcornoral>0
        midc = floor(difcornoral/2);
        CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
    else
        midc = floor(-difcornoral/2);
        CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
    end
    image(CornoralDATnew,'parent',Hshow.Cornoral_axes,'CDataMapping','scaled');axis(Hshow.Cornoral_axes,'off');colormap(Hshow.Cornoral_axes,'gray')
    set(Hshow.Cornoral_axes,'Clim',[30 85]);
else    
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(1)-colmaxt(2));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        CornoralDAT = DOUTSHOWcornoral(:,:,Corval);
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        CornoralDAT = DOUTSHOWcornoralPOS(:,:,Corval);
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS        
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        CornoralDAT = DOUTSHOWcornoralNEG(:,:,Corval);
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG        
    end
    sizeCornoral = size(CornoralDAT);
    difcornoral = sizeCornoral(1)-sizeCornoral(2);
    CornoralDATnew = ones(max(sizeCornoral))*-dshownum;
    if difcornoral>0
        midc = floor(difcornoral/2);
        CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
    else
        midc = floor(-difcornoral/2);
        CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
    end
    image(CornoralDATnew,'parent',Hshow.Cornoral_axes,'CDataMapping','scaled');axis(Hshow.Cornoral_axes,'off');
    if Showtype1
        set(Hshow.Cornoral_axes,'Clim',[-dshownum*2 dshownum*2]);
        colormap(COLMAPBOTH)
    elseif Showtype2
        set(Hshow.Cornoral_axes,'Clim',[-dshownum*1 dshownum*2]);
        colormap(COLMAPPOS)        
    elseif Showtype3
        set(Hshow.Cornoral_axes,'Clim',[-dshownum*2 dshownum*1]);
        colormap(COLMAPNEG)        
    end
end
end
function ASmapshowAxiConSlid(varargin)
Hshow = varargin{3};
inputdir = get(Hshow.IO_inputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
Axival = get(Hshow.Axial_contslide,'val');
Axival = uint8(Axival);
corpos = [1.0,1.0,double(Axival)];
mnipos = cor2mni(corpos,Hshow.Vbg.mat);
MNI = mnipos(3);
set(Hshow.Axial_coredit,'string',num2str(Axival));
set(Hshow.Axial_mniedit,'string',num2str(MNI));
set(Hshow.Axial_cortext,'string',['cor:z=',num2str(Axival)]);
set(Hshow.Axial_mnitext,'string',['mni:z=',num2str(MNI)]);

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');

if isempty(inputdir)
    load([Hshow.pat,filesep,'TempForOrigShow',filesep,'Z',num2str(Axival),'.mat']);
    AxialDAT = rot90(Slice);
    sizeaxial = size(AxialDAT);
    difaxis = sizeaxial(1)-sizeaxial(2);
    AxialDATnew = ones(max(sizeaxial))*30;
    if difaxis>0
        midc = floor(difaxis/2);
        AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
    else
        midc = floor(-difaxis/2);
        AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
    end
    image(AxialDATnew,'parent',Hshow.Axial_axes,'CDataMapping','scaled');axis(Hshow.Axial_axes,'off');colormap(Hshow.Axial_axes,'gray')
    set(Hshow.Axial_axes,'Clim',[30 85]);
else
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(1)-colmaxt(2));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        AxialDAT = DOUTSHOWaxial(:,:,Axival);
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        AxialDAT = DOUTSHOWaxialPOS(:,:,Axival);
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS        
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        AxialDAT = DOUTSHOWaxialNEG(:,:,Axival);
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG        
    end
    sizeaxial = size(AxialDAT);
    difaxis = sizeaxial(1)-sizeaxial(2);
    AxialDATnew = ones(max(sizeaxial))*-dshownum;
    if difaxis>0
        midc = floor(difaxis/2);
        AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
    else
        midc = floor(-difaxis/2);
        AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
    end
    image(AxialDATnew,'parent',Hshow.Axial_axes,'CDataMapping','scaled');axis(Hshow.Axial_axes,'off');
    if Showtype1
        set(Hshow.Axial_axes,'Clim',[-dshownum*2 dshownum*2]);
        colormap(COLMAPBOTH)
    elseif Showtype2
        set(Hshow.Axial_axes,'Clim',[-dshownum*1 dshownum*2]);
        colormap(COLMAPPOS)        
    elseif Showtype3
        set(Hshow.Axial_axes,'Clim',[-dshownum*2 dshownum*1]);
        colormap(COLMAPNEG)        
    end
end
end
function ASmapshowSagChange(varargin)
Hshow = varargin{3};
queans = questdlg('mni or cor?','mni or cor?','mni','cor','mni');
inputdir = get(Hshow.IO_inputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
if strcmp(queans,'mni')
    Axivala = get(Hshow.Sagittal_mniedit,'string');
    mnipos = [ str2num(Axivala) 1.0 1.0];
    corpos = mni2cor(mnipos,Hshow.Vbg.mat);
    SagVal = corpos(1);
    MNI = mnipos(1);
    set(Hshow.Sagittal_coredit,'string',num2str(SagVal));
    set(Hshow.Sagittal_cortext,'string',['cor:z=',num2str(SagVal)]);
    set(Hshow.Sagittal_mnitext,'string',['mni:z=',num2str(MNI)]);
else
    Axivala = get(Hshow.Sagittal_coredit,'string');
    corpos = [str2num(Axivala) 1.0 1.0 ];
    mnipos = cor2mni(corpos,Hshow.Vbg.mat);
    SagVal = corpos(1);
    MNI = mnipos(1);
    set(Hshow.Sagittal_mniedit,'string',num2str(MNI));
    set(Hshow.Sagittal_cortext,'string',['cor:z=',num2str(SagVal)]);
    set(Hshow.Sagittal_mnitext,'string',['mni:z=',num2str(MNI)]);
end
set(Hshow.Sagittal_contslide,'val',SagVal);
% inputdir = get(Hshow.IO_inputed,'string');
% outdir = get(Hshow.IO_Outputed,'string');
% SagVal = get(Hshow.Sagittal_contslide,'val');
% SagVal = uint8(SagVal);
% corpos = [double(SagVal),1.0,1.0];
% mnipos = cor2mni(corpos,Hshow.Vbg.mat);
% MNI = mnipos(1);
% set(Hshow.Sagittal_coredit,'string',num2str(SagVal));
% set(Hshow.Sagittal_mniedit,'string',num2str(MNI));
% set(Hshow.Sagittal_cortext,'string',['cor:x=',num2str(SagVal)]);
% set(Hshow.Sagittal_mnitext,'string',['mni:x=',num2str(MNI)]);

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');

if isempty(inputdir)
    load([Hshow.pat,filesep,'TempForOrigShow',filesep,'X',num2str(SagVal),'.mat']);
    SagittalDAT = rot90(Slice);
    sizeSagittal = size(SagittalDAT);
    difsagittal = sizeSagittal(1)-sizeSagittal(2);
    SagittalDATnew = ones(max(sizeSagittal))*30;
    if difsagittal>0
        midc = floor(difsagittal/2);
        SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
    else
        midc = floor(-difsagittal/2);
        SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
    end
    image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');colormap(Hshow.Sagittal_axes,'gray')
    set(Hshow.Sagittal_axes,'Clim',[30 85]);
else
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(1)-colmaxt(2));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        SagittalDAT = DOUTSHOWsagittal(:,:,SagVal);
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        SagittalDAT = DOUTSHOWsagittalPOS(:,:,SagVal);
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        SagittalDAT = DOUTSHOWsagittalNEG(:,:,SagVal);
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG
    end
    sizeSagittal = size(SagittalDAT);
    difsagittal = sizeSagittal(1)-sizeSagittal(2);
    SagittalDATnew = ones(max(sizeSagittal))*-dshownum;
    if difsagittal>0
        midc = floor(difsagittal/2);
        SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
    else
        midc = floor(-difsagittal/2);
        SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
    end
    image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');
    if Showtype1
        set(Hshow.Sagittal_axes,'Clim',[-dshownum*2 dshownum*2]);
        colormap(COLMAPBOTH)
    elseif Showtype2
        set(Hshow.Sagittal_axes,'Clim',[-dshownum*1 dshownum*2]);
        colormap(COLMAPPOS)        
    elseif Showtype3
        set(Hshow.Sagittal_axes,'Clim',[-dshownum*2 dshownum*1]);
        colormap(COLMAPNEG)        
    end
end

end
function ASmapshowCorChange(varargin)
Hshow = varargin{3};
queans = questdlg('mni or cor?','mni or cor?','mni','cor','mni');
inputdir = get(Hshow.IO_inputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
if strcmp(queans,'mni')
    Axivala = get(Hshow.Cornoral_mniedit,'string');
    mnipos = [1.0  str2num(Axivala) 1.0];
    corpos = mni2cor(mnipos,Hshow.Vbg.mat);
    Corval = corpos(2);
    MNI = mnipos(2);
    set(Hshow.Cornoral_coredit,'string',num2str(Corval));
    set(Hshow.Cornoral_cortext,'string',['cor:z=',num2str(Corval)]);
    set(Hshow.Cornoral_mnitext,'string',['mni:z=',num2str(MNI)]);
else
    Axivala = get(Hshow.Cornoral_coredit,'string');
    corpos = [1.0  str2num(Axivala) 1.0];
    mnipos = cor2mni(corpos,Hshow.Vbg.mat);
    Corval = corpos(2);
    MNI = mnipos(2);
    set(Hshow.Cornoral_mniedit,'string',num2str(MNI));
    set(Hshow.Cornoral_cortext,'string',['cor:z=',num2str(Corval)]);
    set(Hshow.Cornoral_mnitext,'string',['mni:z=',num2str(MNI)]);
end
set(Hshow.Cornoral_contslide,'val',Corval);
% Axival = get(Hshow.Axial_contslide,'val');
% Axival = uint8(Axival);
% corpos = [1.0,1.0,double(Axival)];
% mnipos = cor2mni(corpos,Hshow.Vbg.mat);
% MNI = mnipos(3);
% set(Hshow.Axial_coredit,'string',num2str(Axival));
% set(Hshow.Axial_mniedit,'string',num2str(MNI));
% set(Hshow.Axial_cortext,'string',['cor:z=',num2str(Axival)]);
% set(Hshow.Axial_mnitext,'string',['mni:z=',num2str(MNI)]);
% 
% inputdir = get(Hshow.IO_inputed,'string');
% Corval = get(Hshow.Cornoral_contslide,'val');
% Corval = uint8(Corval);
% corpos = [1.0,double(Corval),1.0];
% mnipos = cor2mni(corpos,Hshow.Vbg.mat);
% MNI = mnipos(2);
% set(Hshow.Cornoral_coredit,'string',num2str(Corval));
% set(Hshow.Cornoral_mniedit,'string',num2str(MNI));
% set(Hshow.Cornoral_cortext,'string',['cor:y=',num2str(Corval)]);
% set(Hshow.Cornoral_mnitext,'string',['mni:y=',num2str(MNI)]);

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');

if isempty(inputdir)
    load([Hshow.pat,filesep,'TempForOrigShow',filesep,'Y',num2str(Corval),'.mat']);
    CornoralDAT = rot90(Slice);
    sizeCornoral = size(CornoralDAT);
    difcornoral = sizeCornoral(1)-sizeCornoral(2);
    CornoralDATnew = ones(max(sizeCornoral))*30;
    if difcornoral>0
        midc = floor(difcornoral/2);
        CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
    else
        midc = floor(-difcornoral/2);
        CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
    end
    image(CornoralDATnew,'parent',Hshow.Cornoral_axes,'CDataMapping','scaled');axis(Hshow.Cornoral_axes,'off');colormap(Hshow.Cornoral_axes,'gray')
    set(Hshow.Cornoral_axes,'Clim',[30 85]);
else    
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(1)-colmaxt(2));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        CornoralDAT = DOUTSHOWcornoral(:,:,Corval);
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        CornoralDAT = DOUTSHOWcornoralPOS(:,:,Corval);
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS        
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        CornoralDAT = DOUTSHOWcornoralNEG(:,:,Corval);
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG        
    end
    sizeCornoral = size(CornoralDAT);
    difcornoral = sizeCornoral(1)-sizeCornoral(2);
    CornoralDATnew = ones(max(sizeCornoral))*-dshownum;
    if difcornoral>0
        midc = floor(difcornoral/2);
        CornoralDATnew(:,midc+1:midc+sizeCornoral(2)) = CornoralDAT;
    else
        midc = floor(-difcornoral/2);
        CornoralDATnew(midc+1:midc+sizeCornoral(1),:) = CornoralDAT;
    end
    image(CornoralDATnew,'parent',Hshow.Cornoral_axes,'CDataMapping','scaled');axis(Hshow.Cornoral_axes,'off');
    if Showtype1
        set(Hshow.Sagittal_axes,'Clim',[-dshownum*2 dshownum*2]);
        colormap(COLMAPBOTH)
    elseif Showtype2
        set(Hshow.Sagittal_axes,'Clim',[-dshownum*1 dshownum*2]);
        colormap(COLMAPPOS)        
    elseif Showtype3
        set(Hshow.Sagittal_axes,'Clim',[-dshownum*2 dshownum*1]);
        colormap(COLMAPNEG)        
    end
end

end
function ASmapshowAxiChange(varargin)
Hshow = varargin{3};
queans = questdlg('mni or cor?','mni or cor?','mni','cor','mni');
inputdir = get(Hshow.IO_inputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
if strcmp(queans,'mni')
    Axivala = get(Hshow.Axial_mniedit,'string');
    mnipos = [1.0 1.0 str2num(Axivala)];
    corpos = mni2cor(mnipos,Hshow.Vbg.mat);
    Axival = corpos(3);
    MNI = mnipos(3);
    set(Hshow.Axial_coredit,'string',num2str(Axival));
    set(Hshow.Axial_cortext,'string',['cor:z=',num2str(Axival)]);
    set(Hshow.Axial_mnitext,'string',['mni:z=',num2str(MNI)]);
else
    Axivala = get(Hshow.Axial_coredit,'string');
    corpos = [1.0 1.0 str2num(Axivala)];
    mnipos = cor2mni(corpos,Hshow.Vbg.mat);
    Axival = corpos(3);
    MNI = mnipos(3);
    set(Hshow.Axial_mniedit,'string',num2str(MNI));
    set(Hshow.Axial_cortext,'string',['cor:z=',num2str(Axival)]);
    set(Hshow.Axial_mnitext,'string',['mni:z=',num2str(MNI)]);
end
set(Hshow.Axial_contslide,'val',Axival);
% Axival = get(Hshow.Axial_contslide,'val');
% Axival = uint8(Axival);
% corpos = [1.0,1.0,double(Axival)];
% mnipos = cor2mni(corpos,Hshow.Vbg.mat);
% MNI = mnipos(3);
% set(Hshow.Axial_coredit,'string',num2str(Axival));
% set(Hshow.Axial_mniedit,'string',num2str(MNI));
% set(Hshow.Axial_cortext,'string',['cor:z=',num2str(Axival)]);
% set(Hshow.Axial_mnitext,'string',['mni:z=',num2str(MNI)]);

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');

if isempty(inputdir)
    load([Hshow.pat,filesep,'TempForOrigShow',filesep,'Z',num2str(Axival),'.mat']);
    AxialDAT = rot90(Slice);
    sizeaxial = size(AxialDAT);
    difaxis = sizeaxial(1)-sizeaxial(2);
    AxialDATnew = ones(max(sizeaxial))*30;
    if difaxis>0
        midc = floor(difaxis/2);
        AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
    else
        midc = floor(-difaxis/2);
        AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
    end
    image(AxialDATnew,'parent',Hshow.Axial_axes,'CDataMapping','scaled');axis(Hshow.Axial_axes,'off');colormap(Hshow.Axial_axes,'gray')
    set(Hshow.Axial_axes,'Clim',[30 85]);
else
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(1)-colmaxt(2));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        AxialDAT = DOUTSHOWaxial(:,:,Axival);
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        AxialDAT = DOUTSHOWaxialPOS(:,:,Axival);
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS        
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        AxialDAT = DOUTSHOWaxialNEG(:,:,Axival);
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG        
    end
    sizeaxial = size(AxialDAT);
    difaxis = sizeaxial(1)-sizeaxial(2);
    AxialDATnew = ones(max(sizeaxial))*-dshownum;
    if difaxis>0
        midc = floor(difaxis/2);
        AxialDATnew(:,midc+1:midc+sizeaxial(2)) = AxialDAT;
    else
        midc = floor(-difaxis/2);
        AxialDATnew(midc+1:midc+sizeaxial(1),:) = AxialDAT;
    end
    image(AxialDATnew,'parent',Hshow.Axial_axes,'CDataMapping','scaled');axis(Hshow.Axial_axes,'off');
    if Showtype1
        set(Hshow.Axial_axes,'Clim',[-dshownum*2 dshownum*2]);
        colormap(COLMAPBOTH)
    elseif Showtype2
        set(Hshow.Axial_axes,'Clim',[-dshownum*1 dshownum*2]);
        colormap(COLMAPPOS)        
    elseif Showtype3
        set(Hshow.Axial_axes,'Clim',[-dshownum*2 dshownum*1]);
        colormap(COLMAPNEG)        
    end
end
end
function ASmapshowSagPrintAll(varargin)
Hshow = varargin{3};
inputdir = get(Hshow.IO_inputed,'string');
templatedir = get(Hshow.IO_tempputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');
Pval = get(Hshow.Other_PthE,'string');
Types = get(Hshow.Other_typeE,'string');
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Clusnum = get(Hshow.Other_Cnum,'string');
SliceSag = get(Hshow.Sagittal_contslide,'val');
SliceCor = get(Hshow.Cornoral_contslide,'val');
SliceAxi = get(Hshow.Axial_contslide,'val');
DOF = get(Hshow.Other_dofE,'string');

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);

DOFv = str2num(DOF);
if strcmpi(Types,'F')&&length(DOFv)<2
    error('wrong Degree of Freedom for F stastical analysis');
end
dof1 = DOFv(1);
if length(DOFv)<2
    dof20 = [];
else
    dof20 = DOFv(2);
end
pthrv = str2num(Pval);
colmaxt = [str2num(Thrval1),str2num(Thrval2)];
CLUSNUM = str2num(Clusnum);
% [DOUTSHOWaxial,DOUTSHOWcornoral,DOUTSHOWsagittal,dshownum]=AS_MapFCshow3sliceU(outdir,inputdir,templatedir,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM);
% SagShowmap = DOUTSHOWsagittal(:,:,SliceSag);
% CorShowmap = DOUTSHOWcornoral(:,:,SliceCor);
% AxiShowmap = DOUTSHOWaxial(:,:,SliceAxi);
if strcmpi(Types,'P')
    ShowtypeTry = [Showtype1,Showtype2,Showtype3];
    OutnameLab = 'PrintAllSlice_Sagittal';
    AS_MapFCshow_used_simpleUPX(outdir,inputdir,templatedir,OutnameLab,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM,ShowtypeTry)
    
else
    ShowtypeTry = [Showtype1,Showtype2,Showtype3];
    OutnameLab = 'PrintAllSlice_Sagittal';
    AS_MapFCshow_used_simpleUX(outdir,inputdir,templatedir,OutnameLab,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM,ShowtypeTry)
    % AS_MapFCshow_used_simple(outdir,inputdir,OutnameLab,pthrv,labmarkinput,dof1,dof20,colmaxt,CLUSNUM)
end

end
function ASmapshowCorPrintAll(varargin)
Hshow = varargin{3};
inputdir = get(Hshow.IO_inputed,'string');
templatedir = get(Hshow.IO_tempputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');
Pval = get(Hshow.Other_PthE,'string');
Types = get(Hshow.Other_typeE,'string');
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Clusnum = get(Hshow.Other_Cnum,'string');
SliceSag = get(Hshow.Sagittal_contslide,'val');
SliceCor = get(Hshow.Cornoral_contslide,'val');
SliceAxi = get(Hshow.Axial_contslide,'val');
DOF = get(Hshow.Other_dofE,'string');

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);

DOFv = str2num(DOF);
if strcmpi(Types,'F')&&length(DOFv)<2
    error('wrong Degree of Freedom for F stastical analysis');
end
dof1 = DOFv(1);
if length(DOFv)<2
    dof20 = [];
else
    dof20 = DOFv(2);
end
pthrv = str2num(Pval);
colmaxt = [str2num(Thrval1),str2num(Thrval2)];
CLUSNUM = str2num(Clusnum);
% [DOUTSHOWaxial,DOUTSHOWcornoral,DOUTSHOWsagittal,dshownum]=AS_MapFCshow3sliceU(outdir,inputdir,templatedir,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM);
% SagShowmap = DOUTSHOWsagittal(:,:,SliceSag);
% CorShowmap = DOUTSHOWcornoral(:,:,SliceCor);
% AxiShowmap = DOUTSHOWaxial(:,:,SliceAxi);
if strcmpi(Types,'P')
    ShowtypeTry = [Showtype1,Showtype2,Showtype3];
    OutnameLab = 'PrintAllSlice_Cornoral';
    AS_MapFCshow_used_simpleUPY(outdir,inputdir,templatedir,OutnameLab,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM,ShowtypeTry)
    
else
    ShowtypeTry = [Showtype1,Showtype2,Showtype3];
    OutnameLab = 'PrintAllSlice_Cornoral';
    AS_MapFCshow_used_simpleUY(outdir,inputdir,templatedir,OutnameLab,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM,ShowtypeTry)
    % AS_MapFCshow_used_simple(outdir,inputdir,OutnameLab,pthrv,labmarkinput,dof1,dof20,colmaxt,CLUSNUM)
end

end
function ASmapshowAxiPrintAll(varargin)
Hshow = varargin{3};
inputdir = get(Hshow.IO_inputed,'string');
templatedir = get(Hshow.IO_tempputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
Thrval1 = get(Hshow.Other_Threshold1,'string');
Thrval2 = get(Hshow.Other_Threshold2,'string');
Pval = get(Hshow.Other_PthE,'string');
Types = get(Hshow.Other_typeE,'string');
Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
Clusnum = get(Hshow.Other_Cnum,'string');
SliceSag = get(Hshow.Sagittal_contslide,'val');
SliceCor = get(Hshow.Cornoral_contslide,'val');
SliceAxi = get(Hshow.Axial_contslide,'val');
DOF = get(Hshow.Other_dofE,'string');

load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);

DOFv = str2num(DOF);
if strcmpi(Types,'F')&&length(DOFv)<2
    error('wrong Degree of Freedom for F stastical analysis');
end
dof1 = DOFv(1);
if length(DOFv)<2
    dof20 = [];
else
    dof20 = DOFv(2);
end
pthrv = str2num(Pval);
colmaxt = [str2num(Thrval1),str2num(Thrval2)];
CLUSNUM = str2num(Clusnum);
% [DOUTSHOWaxial,DOUTSHOWcornoral,DOUTSHOWsagittal,dshownum]=AS_MapFCshow3sliceU(outdir,inputdir,templatedir,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM);
% SagShowmap = DOUTSHOWsagittal(:,:,SliceSag);
% CorShowmap = DOUTSHOWcornoral(:,:,SliceCor);
% AxiShowmap = DOUTSHOWaxial(:,:,SliceAxi);
if strcmpi(Types,'P')
    ShowtypeTry = [Showtype1,Showtype2,Showtype3];
    OutnameLab = 'PrintAllSlice_Axial';
    AS_MapFCshow_used_simpleUPZ(outdir,inputdir,templatedir,OutnameLab,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM,ShowtypeTry)
    
else
    ShowtypeTry = [Showtype1,Showtype2,Showtype3];
    OutnameLab = 'PrintAllSlice_Axial';
    AS_MapFCshow_used_simpleUZ(outdir,inputdir,templatedir,OutnameLab,pthrv,Types,dof1,dof20,colmaxt,CLUSNUM,ShowtypeTry)
    % AS_MapFCshow_used_simple(outdir,inputdir,OutnameLab,pthrv,labmarkinput,dof1,dof20,colmaxt,CLUSNUM)
end

end
function ASmapshowSagPrintCurrent(varargin)
Hshow = varargin{3};
inputdir = get(Hshow.IO_inputed,'string');
outdir = get(Hshow.IO_Outputed,'string');
SagVal = get(Hshow.Sagittal_contslide,'val');
SagVal = uint8(SagVal);
MNI = mnipos(1);
if isempty(inputdir)
    load([Hshow.pat,filesep,'TempForOrigShow',filesep,'X',num2str(SagVal),'.mat']);
    SagittalDAT = rot90(Slice);
    sizeSagittal = size(SagittalDAT);
    difsagittal = sizeSagittal(1)-sizeSagittal(2);
    SagittalDATnew = ones(max(sizeSagittal))*30;
    if difsagittal>0
        midc = floor(difsagittal/2);
        SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
    else
        midc = floor(-difsagittal/2);
        SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
    end
    
    image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');colormap(Hshow.Sagittal_axes,'gray')
    set(Hshow.Sagittal_axes,'Clim',[30 85]);
else
    load([outdir,filesep,'Showinfo.mat']);
%     DOUTSHOWsagittal
    SagittalDAT = DOUTSHOWsagittal(:,:,SagVal);
    clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    sizeSagittal = size(SagittalDAT);
    difsagittal = sizeSagittal(1)-sizeSagittal(2);
    SagittalDATnew = ones(max(sizeSagittal))*30;
    if difsagittal>0
        midc = floor(difsagittal/2);
        SagittalDATnew(:,midc+1:midc+sizeSagittal(2)) = SagittalDAT;
    else
        midc = floor(-difsagittal/2);
        SagittalDATnew(midc+1:midc+sizeSagittal(1),:) = SagittalDAT;
    end
    
    image(SagittalDATnew,'parent',Hshow.Sagittal_axes,'CDataMapping','scaled');axis(Hshow.Sagittal_axes,'off');colormap(Hshow.Sagittal_axes,'gray')
    set(Hshow.Sagittal_axes,'Clim',[30 85]);
end
end
function ASmapshowCorPrintCurrent(varargin)
Hshow = varargin{3};

end
function ASmapshowAxiPrintCurrent(varargin)
Hshow = varargin{3};

end
function ASmapshowSagPrintRange(varargin)
Hshow = varargin{3};
inputrag = get(Hshow.Sagittal_edrag,'string');
Ranges = str2num(inputrag);
load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
inputdir = get(Hshow.IO_inputed,'string');
SagVal = get(Hshow.Sagittal_contslide,'val');
SagVal = uint8(SagVal);
corpos = [double(SagVal),1.0,1.0];
mnipos = cor2mni(corpos,Hshow.Vbg.mat);
outdir = get(Hshow.IO_Outputed,'string');
% save currenttemp
if isempty(inputdir)
    AS_MapFCshow_multisliceUBG(Hshow.pat,'X',Ranges,outdir,'Multislice',colgray,[30 85],1)
else
    Thrval1 = get(Hshow.Other_Threshold1,'string');
    Thrval2 = get(Hshow.Other_Threshold2,'string');
    Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
    Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
    Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(2)-colmaxt(1));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        DATASHOW = DOUTSHOWsagittal;
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        DATASHOW = DOUTSHOWsagittalPOS;
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS        
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        DATASHOW = DOUTSHOWsagittalNEG;
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG        
    end
    if Showtype1==1
        colormapshow = COLMAPBOTH;
        climsuse = [-dshownum*2 dshownum*2];
    elseif Showtype2==1
        colormapshow = COLMAPPOS;
        climsuse = [-dshownum*1 dshownum*2];        
    elseif Showtype3==1
        colormapshow = COLMAPNEG;
        climsuse = [-dshownum*2 dshownum*1];        
    end
    AS_MapFCshow_multisliceU(DATASHOW,Ranges,outdir,'Multislice',colormapshow,climsuse,1);
end
end
function ASmapshowCorPrintRange(varargin)
Hshow = varargin{3};
inputrag = get(Hshow.Cornoral_edrag,'string');
Ranges = str2num(inputrag);
load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
inputdir = get(Hshow.IO_inputed,'string');
Corval = get(Hshow.Cornoral_contslide,'val');
Corval = uint8(Corval);
corpos = [1.0,double(Corval),1.0];
mnipos = cor2mni(corpos,Hshow.Vbg.mat);
MNI = mnipos(2);
outdir = get(Hshow.IO_Outputed,'string');
if isempty(inputdir)
    AS_MapFCshow_multisliceUBG(Hshow.pat,'Y',Ranges,outdir,'Multislice',colgray,[30 85],1)
else
    Thrval1 = get(Hshow.Other_Threshold1,'string');
    Thrval2 = get(Hshow.Other_Threshold2,'string');
    Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
    Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
    Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(2)-colmaxt(1));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        DATASHOW = DOUTSHOWcornoral;
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        DATASHOW = DOUTSHOWcornoralPOS;
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS        
    elseif Showtype3
        load([outdir,filesep,'Showinfo.mat']);
        DATASHOW = DOUTSHOWcornoralNEG;
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG         
    end
    if Showtype1==1
        colormapshow = COLMAPBOTH;
        climsuse = [-dshownum*2 dshownum*2];
    elseif Showtype2==1
        colormapshow = COLMAPPOS;
        climsuse = [-dshownum*1 dshownum*2];        
    elseif Showtype3==1
        colormapshow = COLMAPNEG;
        climsuse = [-dshownum*2 dshownum*1];        
    end
    AS_MapFCshow_multisliceU(DATASHOW,Ranges,outdir,'Multislice',colormapshow,climsuse,1);
    
%     AS_MapFCshow_multisliceU(DATASHOW,sliceorder,outdir,OutnameLab,colormapshow,climsuse,printlab)
end
end
function ASmapshowAxiPrintRange(varargin)
Hshow = varargin{3};
inputrag = get(Hshow.Axial_edrag,'string');
Ranges = str2num(inputrag);
load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
inputdir = get(Hshow.IO_inputed,'string');
Axival = get(Hshow.Axial_contslide,'val');
Axival = uint8(Axival);
corpos = [1.0,1.0,double(Axival)];
mnipos = cor2mni(corpos,Hshow.Vbg.mat);
outdir = get(Hshow.IO_Outputed,'string');
if isempty(inputdir)
    AS_MapFCshow_multisliceUBG(Hshow.pat,'Z',Ranges,outdir,'Multislice',colgray,[30 85],1)
else
    Thrval1 = get(Hshow.Other_Threshold1,'string');
    Thrval2 = get(Hshow.Other_Threshold2,'string');
    Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
    Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
    Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(2)-colmaxt(1));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        DATASHOW = DOUTSHOWaxial;
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        DATASHOW = DOUTSHOWaxialPOS;
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        DATASHOW = DOUTSHOWaxialNEG;
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG
    end
    if Showtype1==1
        colormapshow = COLMAPBOTH;
        climsuse = [-dshownum*2 dshownum*2];
    elseif Showtype2==1
        colormapshow = COLMAPPOS;
        climsuse = [-dshownum*1 dshownum*2];        
    elseif Showtype3==1
        colormapshow = COLMAPNEG;
        climsuse = [-dshownum*2 dshownum*1];        
    end
    AS_MapFCshow_multisliceU(DATASHOW,Ranges,outdir,'Multislice',colormapshow,climsuse,1);
end
end
function ASmapshowSagRangeshow(varargin)
Hshow = varargin{3};
inputrag = get(Hshow.Sagittal_edrag,'string');
Ranges = str2num(inputrag);
load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);

inputdir = get(Hshow.IO_inputed,'string');
SagVal = get(Hshow.Sagittal_contslide,'val');
SagVal = uint8(SagVal);
corpos = [double(SagVal),1.0,1.0];
mnipos = cor2mni(corpos,Hshow.Vbg.mat);
outdir = get(Hshow.IO_Outputed,'string');
% save currenttemp
if isempty(inputdir)
    AS_MapFCshow_multisliceUBG(Hshow.pat,'X',Ranges,outdir,'Multislice',colgray,[30 85],0)
else
    Thrval1 = get(Hshow.Other_Threshold1,'string');
    Thrval2 = get(Hshow.Other_Threshold2,'string');
    Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
    Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
    Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(2)-colmaxt(1));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        DATASHOW = DOUTSHOWsagittal;
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        DATASHOW = DOUTSHOWsagittalPOS;
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        DATASHOW = DOUTSHOWsagittalNEG;
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG
    end
    if Showtype1==1
        colormapshow = COLMAPBOTH;
        climsuse = [-dshownum*2 dshownum*2];
    elseif Showtype2==1
        colormapshow = COLMAPPOS;
        climsuse = [-dshownum*1 dshownum*2];        
    elseif Showtype3==1
        colormapshow = COLMAPNEG;
        climsuse = [-dshownum*2 dshownum*1];        
    end
    AS_MapFCshow_multisliceU(DATASHOW,Ranges,outdir,'Multislice',colormapshow,climsuse,0);
end
end
function ASmapshowCorRangeshow(varargin)
Hshow = varargin{3};
inputrag = get(Hshow.Cornoral_edrag,'string');
Ranges = str2num(inputrag);
load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);

inputdir = get(Hshow.IO_inputed,'string');
Corval = get(Hshow.Cornoral_contslide,'val');
Corval = uint8(Corval);
corpos = [1.0,double(Corval),1.0];
mnipos = cor2mni(corpos,Hshow.Vbg.mat);
MNI = mnipos(2);
outdir = get(Hshow.IO_Outputed,'string');
if isempty(inputdir)
    AS_MapFCshow_multisliceUBG(Hshow.pat,'Y',Ranges,outdir,'Multislice',colgray,[30 85],0)
else
    Thrval1 = get(Hshow.Other_Threshold1,'string');
    Thrval2 = get(Hshow.Other_Threshold2,'string');
    Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
    Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
    Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(2)-colmaxt(1));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        DATASHOW = DOUTSHOWcornoral;
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        DATASHOW = DOUTSHOWcornoralPOS;
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        DATASHOW = DOUTSHOWcornoralPOS;
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralPOS DOUTSHOWaxialNEG
    end
    if Showtype1==1
        colormapshow = COLMAPBOTH;
        climsuse = [-dshownum*2 dshownum*2];
    elseif Showtype2==1
        colormapshow = COLMAPPOS;
        climsuse = [-dshownum*1 dshownum*2];        
    elseif Showtype3==1
        colormapshow = COLMAPNEG;
        climsuse = [-dshownum*2 dshownum*1];        
    end
    AS_MapFCshow_multisliceU(DATASHOW,Ranges,outdir,'Multislice',colormapshow,climsuse,0);
    
%     AS_MapFCshow_multisliceU(DATASHOW,sliceorder,outdir,OutnameLab,colormapshow,climsuse,printlab)
end
end
function ASmapshowAxiRangeshow(varargin)
Hshow = varargin{3};
inputrag = get(Hshow.Axial_edrag,'string');
Ranges = str2num(inputrag);
load(fullfile(Hshow.pat,'graycolmap.mat'));
colormapshow0 = AFNICOLORMAP(64);
colormapshow(1:32,:) = colormapshow0(1:32,:);
colormapshow(33:96,:) = colgray;
colormapshow(97:128,:) = colormapshow0(33:64,:);
COLMAPBOTH = colormapshow;
COLMAPNEG = colormapshow(1:96,:);
COLMAPPOS = colormapshow(33:128,:);
inputdir = get(Hshow.IO_inputed,'string');
Axival = get(Hshow.Axial_contslide,'val');
Axival = uint8(Axival);
corpos = [1.0,1.0,double(Axival)];
mnipos = cor2mni(corpos,Hshow.Vbg.mat);
outdir = get(Hshow.IO_Outputed,'string');
if isempty(inputdir)
    AS_MapFCshow_multisliceUBG(Hshow.pat,'Z',Ranges,outdir,'Multislice',colgray,[30 85],0)
else
    Thrval1 = get(Hshow.Other_Threshold1,'string');
    Thrval2 = get(Hshow.Other_Threshold2,'string');
    Showtype1 = get(Hshow.Other_ShowTYPE1,'val');
    Showtype2 = get(Hshow.Other_ShowTYPE2,'val');
    Showtype3 = get(Hshow.Other_ShowTYPE3,'val');
    colmaxt = [str2num(Thrval1),str2num(Thrval2)];
    dshownum = abs(colmaxt(2)-colmaxt(1));
    if Showtype1
        load([outdir,filesep,'Showinfo.mat']);
        DATASHOW = DOUTSHOWaxial;
        clear DOUTSHOWsagittal DOUTSHOWcornoral DOUTSHOWaxial
    elseif Showtype2
        load([outdir,filesep,'ShowinfoPOS.mat']);
        DATASHOW = DOUTSHOWaxialPOS;
        clear DOUTSHOWsagittalPOS DOUTSHOWcornoralPOS DOUTSHOWaxialPOS
    elseif Showtype3
        load([outdir,filesep,'ShowinfoNEG.mat']);
        DATASHOW = DOUTSHOWaxialNEG;
        clear DOUTSHOWsagittalNEG DOUTSHOWcornoralNEG DOUTSHOWaxialNEG
    end
    if Showtype1==1
        colormapshow = COLMAPBOTH;
        climsuse = [-dshownum*2 dshownum*2];
    elseif Showtype2==1
        colormapshow = COLMAPPOS;
        climsuse = [-dshownum*1 dshownum*2];        
    elseif Showtype3==1
        colormapshow = COLMAPNEG;
        climsuse = [-dshownum*2 dshownum*1];        
    end
    AS_MapFCshow_multisliceU(DATASHOW,Ranges,outdir,'Multislice',colormapshow,climsuse,0);
end
end
function ASmapshowSagRag_ed(varargin)
Hshow = varargin{3};
inputrag = get(Hshow.Sagittal_edrag,'string');
Ranges = str2num(inputrag);
set(Hshow.Sagittal_edrag,'string',num2str(Ranges));
end
function ASmapshowCorRag_ed(varargin)
Hshow = varargin{3};
inputrag = get(Hshow.Cornoral_edrag,'string');
Ranges = str2num(inputrag);
set(Hshow.Cornoral_edrag,'string',num2str(Ranges));
end
function ASmapshowAxiRag_ed(varargin)
Hshow = varargin{3};
inputrag = get(Hshow.Axial_edrag,'string');
Ranges = str2num(inputrag);
set(Hshow.Axial_edrag,'string',num2str(Ranges));
end