function BCCT_makefinalMask
clc;clear all;close all;
ASMASK.fig = figure('Name','Make final mask for calculation',...       
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'color',[0.95 0.95 0.95],...
    'position',[0.20 0.15 0.5 0.5]);
movegui(ASMASK.fig,'center'); 

ASMASK.OutputText = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.05,0.9,0.1,0.08],...
    'style','text',...
    'string',{'Output','Directory'},...
    'Fontsize',10);
ASMASK.OutputEDIT = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.15,0.9,0.7,0.08],...
    'style','edit',...
    'string','NULL');
ASMASK.OutputSel = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.86,0.9,0.09,0.08],...
    'style','pushbutton',...
    'string','...');


ASMASK.InputTEXTROI = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.05,0.8,0.2,0.08],...
    'style','text',...
    'string','Input ROI directory',...
    'Fontsize',10);

ASMASK.InputSelROI = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.25,0.8,0.2,0.08],...
    'style','pushbutton',...
    'string','SELECT');
ASMASK.ShowROI = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.05,0.3,0.4,0.5],...
    'style','Listbox',...
    'string','...');

ASMASK.InputTEXTDAT = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.55,0.8,0.2,0.08],...
    'style','text',...
    'string','Input DAT',...
    'Fontsize',10);
ASMASK.InputSelDAT = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.75,0.8,0.2,0.08],...
    'style','pushbutton',...
    'string','SELECT');
ASMASK.ShowDAT = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.55,0.3,0.4,0.5],...
    'style','Listbox',...
    'string','...');
ASMASK.Reset = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.2,0.15,0.1,0.1],...
    'style','pushbutton',...
    'string','RESET');

ASMASK.Return = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.45,0.15,0.1,0.1],...
    'style','pushbutton',...
    'string','Return');

ASMASK.Done = uicontrol('parent',ASMASK.fig,...
    'units','normalized',...
    'pos',[0.7,0.15,0.1,0.1],...
    'style','pushbutton',...
    'string','DONE');

set(ASMASK.InputSelDAT,'callback',{@ASMASKSFInDat,ASMASK});
set(ASMASK.InputSelROI,'callback',{@ASMASKSFInROI,ASMASK});
set(ASMASK.OutputSel,'callback',{@ASMASKSFOutSel,ASMASK});
set(ASMASK.Reset,'callback',{@ASMASKSFReset,ASMASK});
set(ASMASK.Done,'callback',{@ASMASKSFDone,ASMASK});
set(ASMASK.Return,'callback',{@ASMASKSFReturn,ASMASK});
end

function ASMASKSFInDat(varargin)
ASMASK = varargin{3};
inputnum = inputdlg('No. of Groups','No. of Groups',1,{'2'});
Inputnum = str2num(inputnum{1});
for i = 1:Inputnum
    dirs = uigetdir(pwd,['Group ',num2str(i)]);
    DirOut{i,1} = dirs;
end
set(ASMASK.ShowDAT,'string',DirOut);
end

function ASMASKSFInROI(varargin)
ASMASK = varargin{3};
[Fname,Patname,filtind] = uigetfile({'*.nii';'*.img';'*.*'},'ROI files','MultiSelect','on');
if iscell(Fname)
    for i = 1:length(Fname)
        DirROI{i,1} = fullfile(Patname,Fname{i});
    end
else
    DirROI{1,1} = fullfile(Patname,Fname);
end
set(ASMASK.ShowROI,'string',DirROI);
end

function ASMASKSFOutSel(varargin)
ASMASK = varargin{3};
dirout = uigetdir(pwd,'Outputdir');
set(ASMASK.OutputEDIT,'string',dirout);
end

function ASMASKSFReset(varargin)
ASMASK = varargin{3};
close(ASMASK.fig);
clear ASMASK;
BCCT_makefinalMask;
end

function ASMASKSFReturn(varargin)
ASMASK = varargin{3};
close(ASMASK.fig);
BCCT_WTA_GUI;
end

function ASMASKSFDone(varargin)
ASMASK = varargin{3};
Outdir = get(ASMASK.OutputEDIT,'string');
Indat = get(ASMASK.ShowDAT,'string');
Inroi = get(ASMASK.ShowROI,'string');
CalPara.Outdir = Outdir;
CalPara.Indat = Indat;
CalPara.Inroi = Inroi;
save(fullfile(Outdir,'CalPara.mat'),'CalPara');
BCCT_makefinalMask_main(CalPara);
close(ASMASK.fig);
clear ASMASK;
uiwait(msgbox('Finished'));
BCCT
end