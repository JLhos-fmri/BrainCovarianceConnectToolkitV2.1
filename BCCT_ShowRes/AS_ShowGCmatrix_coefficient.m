function AS_ShowGCmatrix_coefficient(Ptype)
GC_RES.Ptype = Ptype;
Hsize = get(0,'screensize');
MIDPOINT = [Hsize(3)/2,Hsize(4)/2];
Asize = [100*3,100+40];
MaxSIZE = [Hsize(3) Hsize(4)]*0.8;
factor = MaxSIZE./Asize;
factornew = min(factor);
POSSIZE = Asize*factornew;
GC_RES.fig = figure('position',[MIDPOINT(1)-POSSIZE(1)/2,MIDPOINT(2)-POSSIZE(2)/2,POSSIZE(1)/1,POSSIZE(2)/1],'name','GC Matrix(Using coefficient)');
GC_RES.IO = uibuttongroup('parent',GC_RES.fig,...
    'units','norm',...
    'pos',[0.1,0.55,0.8,0.4]);
GC_RES.threshold = uibuttongroup('parent',GC_RES.fig,...
    'units','norm',...
    'pos',[0.1,0.225,0.8,0.3]);
GC_RES.show = uicontrol('parent',GC_RES.fig,...
    'units','norm',...
    'pos',[0.1,0.05,0.3,0.15],...
    'style','pushbutton',...
    'string','Show');
GC_RES.exit = uicontrol('parent',GC_RES.fig,...
    'units','norm',...
    'pos',[0.6,0.05,0.3,0.15],...
    'style','pushbutton',...
    'string','exit');
%%
GC_RES.outtxt = uicontrol('parent',GC_RES.IO,...
    'units','norm',...
    'pos',[0.05,0.85,0.1,0.1],...
    'style','text',...
    'string','Outputdir');
GC_RES.outedit = uicontrol('parent',GC_RES.IO,...
    'units','norm',...
    'pos',[0.2,0.85,0.6,0.1],...
    'style','edit',...
    'string','NULL');
GC_RES.outsel = uicontrol('parent',GC_RES.IO,...
    'units','norm',...
    'pos',[0.85,0.85,0.1,0.1],...
    'style','pushbutton',...
    'string','...');
GC_RES.IO_inputGC = uicontrol('parent',GC_RES.IO,...
    'units','norm',...
    'pos',[0.1,0.70,0.3,0.1],...
    'style','pushbutton',...
    'string','SelInputMats');
GC_RES.IO_inputGCP = uicontrol('parent',GC_RES.IO,...
    'units','norm',...
    'pos',[0.6,0.70,0.3,0.1],...
    'style','pushbutton',...
    'string','SelPmats');
GC_RES.IO_inputGCdir = uicontrol('parent',GC_RES.IO,...
    'units','norm',...
    'pos',[0.1,0.05,0.3,0.6],...
    'style','listbox');
GC_RES.IO_inputGCPdir = uicontrol('parent',GC_RES.IO,...
    'units','norm',...
    'pos',[0.6,0.05,0.3,0.6],...
    'style','listbox');
% if Ptype==1
    set(GC_RES.IO_inputGCP,'enable','off');
    set(GC_RES.IO_inputGCPdir,'enable','off');
% end
%%
GC_RES.TheType1 = uicontrol('parent',GC_RES.threshold,...
    'units','norm',...
    'pos',[0.1,0.85,0.3,0.1],...
    'style','rad',...
    'string','Residual',...
    'val',1);
GC_RES.TheType2 = uicontrol('parent',GC_RES.threshold,...
    'units','norm',...
    'pos',[0.6,0.85,0.3,0.1],...
    'style','rad',...
    'string','Pvalue',...
    'val',0);
GC_RES.TheType1_change = uicontrol('parent',GC_RES.threshold,...
    'units','norm',...
    'pos',[0.1,0.75,0.3,0.08],...
    'style','pushbutton',...
    'string','ChangeValue');
GC_RES.TheType1_changelist = uicontrol('parent',GC_RES.threshold,...
    'units','norm',...
    'pos',[0.1,0.05,0.3,0.7],...
    'style','listbox');
GC_RES.TheType2_Ptxt = uicontrol('parent',GC_RES.threshold,...
    'units','norm',...
    'pos',[0.6,0.75,0.1,0.08],...
    'style','text',...
    'string','Pval = ');
GC_RES.TheType2_Pedit = uicontrol('parent',GC_RES.threshold,...
    'units','norm',...
    'pos',[0.7,0.75,0.2,0.08],...
    'style','edit',...
    'string','0.05');
GC_RES.TheType2_Ptypebtg = uibuttongroup('parent',GC_RES.threshold,...
    'units','norm',...
    'pos',[0.6,0.05,0.3,0.7]);
GC_RES.TheType2_Ptype1 = uicontrol('parent',GC_RES.TheType2_Ptypebtg,...
    'units','norm',...
    'pos',[0.1,0.775,0.8,0.2],...
    'style','rad',...
    'string','Uncorrected',...
    'val',1);
GC_RES.TheType2_Ptype2 = uicontrol('parent',GC_RES.TheType2_Ptypebtg,...
    'units','norm',...
    'pos',[0.1,0.525,0.8,0.2],...
    'style','rad',...
    'string','FPA',...
    'val',0);
GC_RES.TheType2_Ptype3 = uicontrol('parent',GC_RES.TheType2_Ptypebtg,...
    'units','norm',...
    'pos',[0.1,0.275,0.8,0.2],...
    'style','rad',...
    'string','FDR',...
    'val',0);
GC_RES.TheType2_Ptype4 = uicontrol('parent',GC_RES.TheType2_Ptypebtg,...
    'units','norm',...
    'pos',[0.1,0.025,0.8,0.2],...
    'style','rad',...
    'string','Bonf',...
    'val',0);
% if Ptype==1
    set(GC_RES.TheType1,'enable','off');
    set(GC_RES.TheType1_change,'enable','off');    
    set(GC_RES.TheType2,'enable','off');
    set(GC_RES.TheType2_Ptxt,'enable','off');
    set(GC_RES.TheType2_Pedit,'enable','off');
    set(GC_RES.TheType2_Ptype1,'enable','off')
    set(GC_RES.TheType2_Ptype2,'enable','off')
    set(GC_RES.TheType2_Ptype3,'enable','off')
    set(GC_RES.TheType2_Ptype4,'enable','off')
% end
%%
set(GC_RES.outsel,'callback',{@OutSel,GC_RES});
set(GC_RES.IO_inputGC,'callback',{@InputSelGC,GC_RES});
set(GC_RES.IO_inputGCP,'callback',{@InputSelGCP,GC_RES});
set(GC_RES.TheType1,'callback',{@TheType1,GC_RES});
set(GC_RES.TheType1_change,'callback',{@Type1Change,GC_RES});
set(GC_RES.TheType2,'callback',{@TheType2,GC_RES});
set(GC_RES.TheType2_Ptype1,'callback',{@TheType2Ptype1,GC_RES});
set(GC_RES.TheType2_Ptype2,'callback',{@TheType2Ptype2,GC_RES});
set(GC_RES.TheType2_Ptype3,'callback',{@TheType2Ptype3,GC_RES});
set(GC_RES.TheType2_Ptype4,'callback',{@TheType2Ptype4,GC_RES});
set(GC_RES.exit,'callback',{@Exitfun,GC_RES});
set(GC_RES.show,'callback',{@ShowMat,GC_RES});
end
function OutSel(varargin)
GC_RES = varargin{3};
dirout = uigetdir(pwd);
set(GC_RES.outedit,'string',dirout);
cd(dirout);
end
function InputSelGC(varargin)
GC_RES = varargin{3};
[nam,path,ext] = uigetfile('*.mat','Matrix of Connection','Multiselect','on');
if iscell(nam)
    for i = 1:length(nam)
        listbox{i,1} = fullfile(path,nam{i});
    end
else
    listbox{1,1} = fullfile(path,nam);
end
set(GC_RES.IO_inputGCdir,'string',listbox);
if GC_RES.Ptype~=1
    set(GC_RES.IO_inputGCP,'enable','on')
    set(GC_RES.IO_inputGCPdir,'enable','on')
end
for i = 1:size(listbox,1)
    TempV = load(listbox{i,1});
    ListVarNameT = fieldnames(TempV);
    ListVarName{i,1} = ListVarNameT{1};
    matrixT = getfield(TempV,ListVarNameT{1});
    matrixT(isnan(matrixT)) = 0;
    matrixT(isinf(matrixT)) = 0;
    maxV(i,1) = max(matrixT(:));
    minV(i,1) = min(matrixT(:));
    clear TempV ListVarNameT matrxiT
end
for i = 1:size(listbox,1)
    listboxthr{i,1} = [ListVarName{i,1},',max:',num2str(maxV(i)),';min:',num2str(minV(i))];
end
set(GC_RES.TheType1_changelist,'string',listboxthr);
%%
set(GC_RES.TheType1,'enable','on');
set(GC_RES.TheType1_change,'enable','on');
if GC_RES.Ptype~=1
    set(GC_RES.TheType2,'enable','on');
    set(GC_RES.TheType2_Ptxt,'enable','on');
    set(GC_RES.TheType2_Pedit,'enable','on');
%     set(GC_RES.TheType2_Ptype1,'enable','on')
%     set(GC_RES.TheType2_Ptype2,'enable','on')
%     set(GC_RES.TheType2_Ptype3,'enable','on')
%     set(GC_RES.TheType2_Ptype4,'enable','on')
end
%%

end
function InputSelGCP(varargin)
GC_RES = varargin{3};
listbox1 = get(GC_RES.IO_inputGCdir,'string');
for i = 1:size(listbox1,1)
    [P1,P2,P3] = fileparts(listbox1{i,1});
    [nam,path,ext] = uigetfile('*.mat',['Pval for ',P2]);
    listboxP{i,1} = fullfile(path,nam);
end
set(GC_RES.IO_inputGCPdir,'string',listboxP);
end
function TheType1(varargin)
GC_RES = varargin{3};
set(GC_RES.TheType1,'val',1);
set(GC_RES.TheType2,'val',0);
set(GC_RES.TheType1_change,'enable','on');
set(GC_RES.TheType1_changelist,'enable','on');
set(GC_RES.TheType2_Ptxt,'enable','off');
set(GC_RES.TheType2_Pedit,'enable','off');
set(GC_RES.TheType2_Ptype1,'enable','off')
set(GC_RES.TheType2_Ptype2,'enable','off')
set(GC_RES.TheType2_Ptype3,'enable','off')
set(GC_RES.TheType2_Ptype4,'enable','off')
end
function Type1Change(varargin)
GC_RES = varargin{3};
% save test111
strings = get(GC_RES.TheType1_changelist,'string');
for i = 1:size(strings,1)
    strtemp = strings{i};
    K = find(strtemp==':');
    k0 = find(strtemp==';');
    maxv = strtemp(K(1)+1:k0-1);
    MAXV = str2num(maxv);
    minv = strtemp(K(2)+1:end);
    MINV = str2num(minv);
    K1 = find(strtemp==',');
    varname = strtemp(1:K1-1);
    Answtemp = inputdlg({['threshold setup: ',varname,',max(pos):',maxv],['maxv(neg):',minv]},['threshold setup: ',varname],1,{num2str(MAXV/2),num2str(MINV/2)});
    LISTSTRING{i,1} = [varname,',Threshold,[',Answtemp{1},',',Answtemp{2},'],max:',maxv];
end
set(GC_RES.TheType1_changelist,'string',LISTSTRING);
end
function TheType2(varargin)
GC_RES = varargin{3};
set(GC_RES.TheType1,'val',0);
set(GC_RES.TheType2,'val',1);
set(GC_RES.TheType1_change,'enable','off');
set(GC_RES.TheType1_changelist,'enable','off');
set(GC_RES.TheType2_Ptxt,'enable','on');
set(GC_RES.TheType2_Pedit,'enable','on');
% set(GC_RES.TheType2_Ptype1,'enable','on')
% set(GC_RES.TheType2_Ptype2,'enable','on')
% set(GC_RES.TheType2_Ptype3,'enable','on')
% set(GC_RES.TheType2_Ptype4,'enable','on')

end
function TheType2Ptype1(varargin)
GC_RES = varargin{3};
set(GC_RES.TheType2Ptype1,'val',1);
set(GC_RES.TheType2Ptype2,'val',0);
set(GC_RES.TheType2Ptype3,'val',0);
set(GC_RES.TheType2Ptype4,'val',0);
end
function TheType2Ptype2(varargin)
GC_RES = varargin{3};
set(GC_RES.TheType2Ptype1,'val',0);
set(GC_RES.TheType2Ptype2,'val',1);
set(GC_RES.TheType2Ptype3,'val',0);
set(GC_RES.TheType2Ptype4,'val',0);
end
function TheType2Ptype3(varargin)
GC_RES = varargin{3};
set(GC_RES.TheType2Ptype1,'val',0);
set(GC_RES.TheType2Ptype2,'val',0);
set(GC_RES.TheType2Ptype3,'val',1);
set(GC_RES.TheType2Ptype4,'val',0);
end
function TheType2Ptype4(varargin)
GC_RES = varargin{3};
set(GC_RES.TheType2Ptype1,'val',0);
set(GC_RES.TheType2Ptype2,'val',0);
set(GC_RES.TheType2Ptype3,'val',0);
set(GC_RES.TheType2Ptype4,'val',1);
end
function Exitfun(varargin)
GC_RES = varargin{3};
close(GC_RES.fig);
end
function ShowMat(varargin)
GC_RES = varargin{3};
Matlist = get(GC_RES.IO_inputGCdir,'string');
for i = 1:size(Matlist,1)
    TempV = load(Matlist{i,1});
    ListVarNameT = fieldnames(TempV);
    ListVarName{i,1} = ListVarNameT{1};
    matrixT = getfield(TempV,ListVarNameT{1});
    matrixT(isnan(matrixT)) = 0;
    matrixT(isinf(matrixT)) = 0;
    MatrixShow{i,1} = matrixT;
    clear TempV ListVarNameT matrxiT
end
save test2
if GC_RES.Ptype==1||(GC_RES.Ptype==2&&get(GC_RES.TheType1,'val')==1) % With No Pval 
    strings = get(GC_RES.TheType1_changelist,'string');
    for i = 1:size(strings,1)
        strtemp = strings{i};
        Ktemp1 = find(strtemp=='[');
        Ktemp2 = find(strtemp==']');
        if isempty(Ktemp1)
            K = find(strtemp==':');
            k0 = find(strtemp==';');
            maxv = strtemp(K(1)+1:k0-1);
            MAXV = str2num(maxv);
            minv = strtemp(K(2)+1:end);
            MINV = str2num(minv);
            thr(i,:) = [MAXV/2,MINV];
        else
            thr(i,:) = str2num(strtemp(Ktemp1:Ktemp2));
        end
        MATSHOW = MatrixShow{i,1};
        MATSHOW(MATSHOW<thr(i,1)&MATSHOW>thr(i,2))=0;
        MATSHOWUse{i,1} = MATSHOW;
    end
else
    PTYPEINSHOW1 = get(GC_RES.TheType2_Ptype1,'val');
    PTYPEINSHOW2 = get(GC_RES.TheType2_Ptype2,'val');
    PTYPEINSHOW3 = get(GC_RES.TheType2_Ptype3,'val');
    PTYPEINSHOW4 = get(GC_RES.TheType2_Ptype4,'val');
    pvalst = get(GC_RES.TheType2_Pedit,'string');
    pvals = str2num(pvalst);
    P_LISTBOX = get(GC_RES.IO_inputGCPdir,'string');
    for i = 1:size(P_LISTBOX,1)        
        TempV = load(P_LISTBOX{i,1});
        P_ListVarNameT = fieldnames(TempV);
        P_ListVarName{i,1} = P_ListVarNameT{1};
        P_matrixT = getfield(TempV,P_ListVarNameT{1});
        PMATRIXSHOW{i,1} = P_matrixT;
    end
    if PTYPEINSHOW1
        for i = 1:size(PMATRIXSHOW,1)
            Ptempmat = PMATRIXSHOW{i,1};
            MATSHOW = MatrixShow{i,1};
            MATSHOW(Ptempmat<(1-pvals)&Ptempmat>pvals)=0;
            MATSHOWUse{i,1} = MATSHOW;
        end
    elseif PTYPEINSHOW2 % FPA comming soon
        uiwait(msgbox('FPA correction, comming soon'));
    elseif PTYPEINSHOW3 % FDR comming soon
        uiwait(msgbox('FDR correction, comming soon'));        
    elseif PTYPEINSHOW4 % Bonf comming soon
        uiwait(msgbox('Bonf correction, comming soon'));        
    end
end
for i = 1:size(MATSHOWUse,1)
    
    Hsize = get(0,'screensize');
    MIDPOINT = [Hsize(3)/2,Hsize(4)/2];
    Asize = [100,100];
    MaxSIZE = [Hsize(3) Hsize(4)]*0.8;
    factor = MaxSIZE./Asize;
    factornew = min(factor);
    POSSIZE = Asize*factornew;
    Hfig{i} = figure('position',[MIDPOINT(1)-POSSIZE(1)/2,MIDPOINT(2)-POSSIZE(2)/2,POSSIZE(1)/1,POSSIZE(2)/1],'name',ListVarName{i,1});
    MATSHOWUsed_S = MATSHOWUse{i,1};
    nsize = size(MATSHOWUsed_S);
    MATSHOWUsed_S(nsize(2),1:nsize(2)) = 0;
    imagesc(MATSHOWUsed_S,[-max(abs(MATSHOWUsed_S(:))),max(abs(MATSHOWUsed_S(:)))]*1.5);title(ListVarName{i,1});
    hold on;
    load('COL_Coef.mat');
    colormap(COL_Coef);
    for j = 1:nsize(2)
        
        plot([0.5,nsize(2)+0.5],[j-0.5,j-0.5],'k--','linewidth',1);
        plot([j-0.5,j-0.5],[0.5,nsize(2)+0.5],'k--','linewidth',1);
        
        plot([j-0.5,nsize(2)+0.5],[j+0.5,j+0.5],'k','linewidth',2);
        plot([j+0.5,j+0.5],[0.5,j+1.5],'k','linewidth',2);
    end
end

end