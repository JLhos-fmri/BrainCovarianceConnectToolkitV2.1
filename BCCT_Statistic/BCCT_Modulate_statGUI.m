function BCCT_Modulate_statGUI
D.fig = figure('Name','BCCT stastic(two groups Modulations)',...            
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'color',[0.95 0.95 0.95],...
    'position',[0.3 0.4 0.4 0.3]);
movegui(D.fig,'center'); 

%%
D.OutputText = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.75,0.1,0.2],...
    'style','text',...
    'string',{'Output','Directory'},...
    'Fontsize',10);
D.OutputEDIT = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.75,0.7,0.2],...
    'style','edit',...
    'string','NULL');
D.OutputSel = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.75,0.09,0.2],...
    'style','pushbutton',...
    'string','...');

%%

D.InputTEXT1 = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.55,0.1,0.2],...
    'style','text',...
    'string',{'Group1 Input','Directory'},...
    'Fontsize',10);
D.InputEDIT1 = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.55,0.7,0.2],...
    'style','edit',...
    'string','NULL');
D.InputSel1 = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.55,0.09,0.2],...
    'style','pushbutton',...
    'string','...');

D.InputTEXT2 = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.05,0.35,0.1,0.2],...
    'style','text',...
    'string',{'Group2 Input','Directory'},...
    'Fontsize',10);
D.InputEDIT2 = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.15,0.35,0.7,0.2],...
    'style','edit',...
    'string','NULL');
D.InputSel2 = uicontrol('parent',D.fig,...
    'units','normalized',...
    'pos',[0.86,0.35,0.09,0.2],...
    'style','pushbutton',...
    'string','...');

%%

D.STType1 = uicontrol('parent',D.fig,...
    'units','norm',...
    'pos',[0.2,0.2,0.1,0.15],...
    'style','text',...
    'string','Interaction');
% D.STType2 = uicontrol('parent',D.fig,...
%     'units','norm',...
%     'pos',[0.5,0.2,0.1,0.15],...
%     'style','rad',...
%     'string','Permutation');
% D.STpermtext = uicontrol('parent',D.fig,...
%     'units','norm',...
%     'pos',[0.65,0.2,0.1,0.15],...
%     'style','text',...
%     'string',{'Permutation', 'times'});
% D.STpermedit = uicontrol('parent',D.fig,...
%     'units','norm',...
%     'pos',[0.75,0.2,0.1,0.15],...
%     'style','edit',...
%     'string','5000');

%%

D.Comp = uicontrol('parent',D.fig,...
    'units','norm',...
    'pos',[0.1,0.1,0.3,0.08],...
    'style','pushbutton',...
    'string','Compute',...
    'Fontsize',10);
D.Exit = uicontrol('parent',D.fig,...
    'units','norm',...
    'pos',[0.6,0.1,0.3,0.08],...
    'style','pushbutton',...
    'string','Exit',...
    'Fontsize',10);


set(D.OutputSel,'callback',{@OutputSel,D});
set(D.InputSel1,'callback',{@InputSel1,D});
set(D.InputSel2,'callback',{@InputSel2,D});

set(D.Exit,'callback',{@Exit,D});
set(D.Comp,'callback',{@calcustat,D});
end

function calcustat(varargin)
D = varargin{3};

Outputdir = get(D.OutputEDIT,'string');
if strcmp(Outputdir,'NULL')
    errordlg('Please Select Outputdir');
    error('Please Select Outputdir');
end
Parameter.Outputdir = Outputdir;
Inputdir1 = get(D.InputEDIT1,'string');
if strcmp(Inputdir1,'NULL')
    errordlg('Please Select Inputdir');
    error('Please Select Inputdir');
end
Parameter.Inputdir1 = Inputdir1;
Inputdir2 = get(D.InputEDIT2,'string');
if strcmp(Inputdir2,'NULL')
    errordlg('Please Select Inputdir');
    error('Please Select Inputdir');
end
Parameter.Inputdir2 = Inputdir2;

% Parameter.Intlab = get(D.STType1,'val');
Parameter.Intlab = 1;
% Parameter.Permlab = get(D.STType2,'val');
% Parameter.PermNum = str2num(get(D.STpermedit,'string'));

outputmat = fullfile(Outputdir,'SetUpparameter.mat');
save(outputmat,'Parameter');
BCCT_Modulation_stat_mainfunc(Parameter)
uiwait(msgbox('Finished'));
end

function Exit(varargin)
D = varargin{3};
close(D.fig);
BCCT;
end


function OutputSel(varargin)
D = varargin{3};
PG = uigetdir(pwd,'Output Directory Selection');
set(D.OutputEDIT,'string',PG)
end
function InputSel1(varargin)
D = varargin{3};
PG = uigetdir(pwd,'Input Directory Selection');
set(D.InputEDIT1,'string',PG)
end
function InputSel2(varargin)
D = varargin{3};
PG = uigetdir(pwd,'Input Directory Selection');
set(D.InputEDIT2,'string',PG)
end