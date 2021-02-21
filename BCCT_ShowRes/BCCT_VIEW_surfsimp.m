function BCCT_VIEW_surfsimp
Hvs.fig = figure('Name','Simple Version for the Surface results viewer',...       
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'color',[0.95 0.95 0.95],...
    'position',[0.3 0.3 0.4 0.4]);
movegui(Hvs.fig,'center'); 

Hvs.BUT1 = uibuttongroup('parent',Hvs.fig,'unit','norm','pos',[0.1,0.8,0.8,0.1]);
Hvs.surftype(1) = uicontrol('parent',Hvs.BUT1,'unit','norm','pos',[0/6,0.1,1/6,0.8],'style','rad','string','fsaverage');
Hvs.surftype(2) = uicontrol('parent',Hvs.BUT1,'unit','norm','pos',[1/6,0.1,1/6,0.8],'style','rad','string','fsaverage3');
Hvs.surftype(3) = uicontrol('parent',Hvs.BUT1,'unit','norm','pos',[2/6,0.1,1/6,0.8],'style','rad','string','fsaverage4');
Hvs.surftype(4) = uicontrol('parent',Hvs.BUT1,'unit','norm','pos',[3/6,0.1,1/6,0.8],'style','rad','string','fsaverage5');
Hvs.surftype(5) = uicontrol('parent',Hvs.BUT1,'unit','norm','pos',[4/6,0.1,1/6,0.8],'style','rad','string','fsaverage6');
Hvs.surftype(6) = uicontrol('parent',Hvs.BUT1,'unit','norm','pos',[5/6,0.1,1/6,0.8],'style','rad','string','fsaverage_sym');

Hvs.output = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.1,0.7,0.2,0.1],'style','text','string','Outputdir');
Hvs.outputed = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.3,0.7,0.5,0.1],'style','edit','string','Null');
Hvs.outputsel = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.8,0.7,0.1,0.1],'style','pushbutton','string','...');



Hvs.dattype(1) = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.1,0.6,0.8,0.1],'style','rad','string','SCN/CaSCN/Modulate for one group or interaction test results');
Hvs.dattype(2) = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.1,0.4,0.8,0.1],'style','rad','string','Permutation Test for Group Comparison');
set(Hvs.dattype(1),'val',1);
set(Hvs.dattype(2),'val',0);
Hvs.stmap = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.1,0.5,0.2,0.1],'style','text','string',{'scn/cascn/modulate','files'});
Hvs.stmaped = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.3,0.5,0.5,0.1],'style','edit','string','Null');
Hvs.stmapsel = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.8,0.5,0.1,0.1],'style','pushbutton','string','...');


Hvs.permmap = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.1,0.3,0.2,0.1],'style','text','string',{'PermTest for ','Group comparison'});
Hvs.permmaped = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.3,0.3,0.5,0.1],'style','edit','string','Null');
Hvs.permmapsel = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.8,0.3,0.1,0.1],'style','pushbutton','string','...');
set(Hvs.permmap,'enable','off');
set(Hvs.permmaped,'enable','off');
set(Hvs.permmapsel,'enable','off');


Hvs.pval = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.1,0.1,0.2,0.1],'style','text','string','Pvalues');
Hvs.pvaled = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.3,0.1,0.2,0.1],'style','edit','string',0.05);

Hvs.show = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.6,0.1,0.1,0.1],'style','pushbutton','string','show');
Hvs.exit = uicontrol('parent',Hvs.fig,'unit','norm','pos',[0.8,0.1,0.1,0.1],'style','pushbutton','string','exit');

set(Hvs.outputsel,'callback',{@Outsel,Hvs});
set(Hvs.dattype(1),'callback',{@dattypesel1,Hvs});
set(Hvs.dattype(2),'callback',{@dattypesel2,Hvs});
set(Hvs.stmapsel,'callback',{@selstatmap,Hvs});
set(Hvs.permmapsel,'callback',{@selpermmap,Hvs});
set(Hvs.show,'callback',{@showres,Hvs});
set(Hvs.exit,'callback',{@exit,Hvs});
end

function Outsel(varargin)
Hvs = varargin{3};
PG = uigetdir(pwd,'Output Directory Selection');
set(Hvs.outputed,'string',PG);
end
function dattypesel1(varargin)
Hvs = varargin{3};
set(Hvs.dattype(1),'val',1);
set(Hvs.dattype(2),'val',0);
set(Hvs.permmap,'enable','off');
set(Hvs.permmaped,'enable','off');
set(Hvs.permmapsel,'enable','off');
set(Hvs.stmap,'enable','on');
set(Hvs.stmaped,'enable','on');
set(Hvs.stmapsel,'enable','on');
end

function dattypesel2(varargin)
Hvs = varargin{3};
set(Hvs.dattype(1),'val',0);
set(Hvs.dattype(2),'val',1);
set(Hvs.permmap,'enable','on');
set(Hvs.permmaped,'enable','on');
set(Hvs.permmapsel,'enable','on');
set(Hvs.stmap,'enable','off');
set(Hvs.stmaped,'enable','off');
set(Hvs.stmapsel,'enable','off');
end

function selstatmap(varargin)
Hvs = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.mgh';'*.txt'},'stat files selection');
set(Hvs.stmaped,'string',fullfile(PathName,FileName));
end
function selpermmap(varargin)
Hvs = varargin{3};
[FileName,PathName,FilterIndex] = uigetfile({'*.mgh';'*.txt'},'perm files selection');
set(Hvs.permmaped,'string',fullfile(PathName,FileName));
end

function showres(varargin)
Hvs = varargin{3};
Outputdir = get(Hvs.outputed,'string');
Parameter.Outputdir = Outputdir;

for i = 1:6
    FStype(i,1) = get(Hvs.surftype(i),'val');
end
Parameter.FSTYPE = find(FStype);
Pval = get(Hvs.pvaled,'string');
Parameter.PVAL = str2num(Pval);
dattype = get(Hvs.dattype(1),'val');
Parameter.DATATYPE = dattype;
% save test
if dattype    
    Parameter.PermTestlab = 0;
    filedir = get(Hvs.stmaped,'string');
    [pat nam ext] = fileparts(filedir);
    lhfilename = filedir;
    rhfilename = fullfile(pat,['r',nam(2:end),'.mgh']);
    Parameter.LHRHfile = {lhfilename,rhfilename};
    
    if strcmp(nam(4:5),'R_') % SCN,mod,interaction
        lhfilename_p = fullfile(pat,['lh.P_',nam(6:end),'.mgh']);
        rhfilename_p = fullfile(pat,['rh.P_',nam(6:end),'.mgh']);
        Parameter.LHRHfile_P = {lhfilename_p,rhfilename_p};
        Parameter.PermPLab = 0;
        Parameter.DATATYPE = 1;
    elseif strcmp(nam(4:5),'T_') % SCN,mod,interaction
        lhfilename_p = fullfile(pat,['lh.P_',nam(6:end),'.mgh']);
        rhfilename_p = fullfile(pat,['rh.P_',nam(6:end),'.mgh']);
        Parameter.LHRHfile_P = {lhfilename_p,rhfilename_p};
        Parameter.PermPLab = 0;
        Parameter.DATATYPE = 1;
    elseif strcmp(nam(4:5),'Z_') % SCN,mod,interaction
        lhfilename_p = fullfile(pat,['lh.P_',nam(6:end),'.mgh']);
        rhfilename_p = fullfile(pat,['rh.P_',nam(6:end),'.mgh']);
        Parameter.LHRHfile_P = {lhfilename_p,rhfilename_p};
        Parameter.PermPLab = 0;
        Parameter.DATATYPE = 1;
    elseif strcmp(nam(4:6),'ROI') % res cascn
        Parameter.DATATYPE = 2;
        if strcmp(nam(13:15),'Net')
            lhfilename_perm = fullfile(pat,['lh.Perm_',nam(4:end),'.mgh']);
            rhfilename_perm = fullfile(pat,['rh.Perm_',nam(4:end),'.mgh']);
            if isempty(dir(lhfilename_perm))
                Parameter.PermPLab = 0;
                uiwait(msgbox('we only represent the unthresholded RESULTS'));
            else
                Parameter.PermPLab = 1;
                Parameter.LHRHfile_Perm = {lhfilename_perm,rhfilename_perm};
            end
        else
            namtotal = nam(4:15);
            lhfilename_p = fullfile(pat,['lh.Pval_',namtotal,'.mgh']);
            rhfilename_p = fullfile(pat,['rh.Pval_',namtotal,'.mgh']);
            lhfilename_perm = fullfile(pat,['lh.Perm_',nam(4:end),'.mgh']);
            rhfilename_perm = fullfile(pat,['rh.Perm_',nam(4:end),'.mgh']);
            Parameter.LHRHfile_P = {lhfilename_p,rhfilename_p};
            if isempty(dir(lhfilename_perm))
                Parameter.PermPLab = 0;
            else
                Parameter.PermPLab = 1;
                Parameter.LHRHfile_Perm = {lhfilename_perm,rhfilename_perm};
            end
        end
    elseif strcmp(nam(4:7),'Coef') % coef cascn
        Parameter.DATATYPE = 3;
        indstart = strfind(nam,'ROI');
        lhfilename_p = fullfile(pat,['lh.Coef_P_',nam(indstart:end),'.mgh']);
        rhfilename_p = fullfile(pat,['rh.Coef_P_',nam(indstart:end),'.mgh']);
        lhfilename_perm = fullfile(pat,['lh.Coef_PermP_',nam(4:end),'.mgh']);
        rhfilename_perm = fullfile(pat,['rh.Coef_PermP_',nam(4:end),'.mgh']);
        Parameter.LHRHfile_P = {lhfilename_p,rhfilename_p};
        if isempty(dir(lhfilename_perm))
            Parameter.PermPLab = 0;
        else
            Parameter.PermPLab = 1;
            Parameter.LHRHfile_Perm = {lhfilename_perm,rhfilename_perm};
        end
    end
else
    filedir = get(Hvs.permmaped,'string');
    [pat nam ext] = fileparts(filedir);
    lhfilename = filedir;
    rhfilename = fullfile(pat,['r',nam(2:end),'.mgh']);
    Parameter.PermTestlab = 1;
    Parameter.PermLHRHfile = {lhfilename,rhfilename};
end
Parameter.name = nam(4:end-4);
save(fullfile(Outputdir,'View_SetUpParameter.mat'),'Parameter');
Surf_Resshow_subfun(Parameter);

end

function exit(varargin)
Hvs = varargin{3};
close(Hvs.fig);
BCCT;
end