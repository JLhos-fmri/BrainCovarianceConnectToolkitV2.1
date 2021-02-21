function BCCT_ShowWTAGUI
has.fig = figure('units','norm','pos',[0.3 0.4 0.4 0.2]);
has.infocol = uicontrol('parent',has.fig,'units','norm','pos',[0.2 0.7 0.2 0.2],'style','pushbutton','string','Collect Information');
has.infocoled = uicontrol('parent',has.fig,'units','norm','pos',[0.4 0.7 0.4 0.2],'style','text','string','Info saved in XXXX');
has.Mapshow = uicontrol('parent',has.fig,'units','norm','pos',[0.2 0.5 0.6 0.2],'style','pushbutton','string','Show Pattern');
has.Radshow = uicontrol('parent',has.fig,'units','norm','pos',[0.2 0.3 0.6 0.2],'style','pushbutton','string','Show Radar chart');
has.Exit = uicontrol('parent',has.fig,'units','norm','pos',[0.2 0.1 0.6 0.2],'style','pushbutton','string','Exit');
set(has.infocol,'callback',{@infocollect,has});
set(has.Mapshow,'callback',{@Mapshow,has});
set(has.Radshow,'callback',{@Radshow,has});
set(has.Exit,'callback',{@Exithas,has});
end
function infocollect(varargin)
has = varargin{3};
% movegui(has.fig,'west');
aswtashowinfo(has);
% movegui(has.fig,'center');
end
function Mapshow(varargin)
has = varargin{3};
% close(has.fig);
INFODIR = get(has.infocoled,'string');
marks = find(INFODIR==':');
if ~isempty(marks)
    INFODIRS = INFODIR(marks+1:end);
else
    return;
end
load(INFODIRS);
[pat,nam,ext] = fileparts(which('BCCT_ShowWTAGUI.m'));
bgpath = fullfile(pat,'mni_icbm152_t1_tal_nlin_asym_09a.nii');
[vbg,dbg] = Dynamic_read_dir_NIFTI(bgpath);
dbgre = reshape(dbg,vbg.dim(1),vbg.dim(2),vbg.dim(3));
Info.bgpath = bgpath;
Hsize = get(0,'screensize');
targetnum = Info.targetnum;
ASSHOWMAP.fig = figure('units','norm','pos',[0.01 0.2 0.1 0.6]);
ASSHOWMAP.Info = Info;
ASSHOWMAP.axi = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.05 0.9 0.3 0.08],'style','rad','string','Axi','val',1);
ASSHOWMAP.cor = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.35 0.9 0.3 0.08],'style','rad','string','Cor','val',0);
ASSHOWMAP.sag = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.65 0.9 0.3 0.08],'style','rad','string','Sag','val',0);
[vtar dtar] = Dynamic_read_dir_NIFTI(Info.target);
dtarind = unique(dtar);
dtarre = reshape(dtar,vtar.dim(1),vtar.dim(2),vtar.dim(3));
indnew = find(dtarre>0);
[indnewx, indnewy, indnewz] = ind2sub(vtar.dim,indnew);
indnewsub = [indnewx, indnewy, indnewz];
indnewsubT = trans2solution(indnewsub,vtar.mat,vbg.mat);
diffx = abs(vtar.mat(1,1)/vbg.mat(1,1))/2*abs(vbg.mat(1,1));
diffy = abs(vtar.mat(2,2)/vbg.mat(2,2))/2*abs(vbg.mat(2,2));
diffz = abs(vtar.mat(3,3)/vbg.mat(3,3))/2*abs(vbg.mat(3,3));
MATOUTT = zeros(vbg.dim);
for i = 1:size(indnewsub,1)
    indu = indnewsubT(i,:);
    induextendx = round(indu(1)-diffx):round(indu(1)+diffx);
    induextendy = round(indu(2)-diffy):round(indu(2)+diffy);
    induextendz = round(indu(3)-diffz):round(indu(3)+diffz);
    MATOUTT(induextendx,induextendy,induextendz) = dtarre(indnew(i));
end
INDSHOW = find(MATOUTT>0);
[ix iy iz] = ind2sub(vbg.dim,INDSHOW);
rangx = unique(ix);
rangy = unique(iy);
rangz = unique(iz);
maxminx = [min(rangx),max(rangx)];
maxminy = [min(rangy),max(rangy)];
maxminz = [min(rangz),max(rangz)];
ASSHOWMAP.maxminx = maxminx;
ASSHOWMAP.maxminy = maxminy;
ASSHOWMAP.maxminz = maxminz;
inmixdir = fullfile(Info.input,'mixedMaxID.nii');
[vmid dmid] = Dynamic_read_dir_NIFTI(inmixdir);
dmidre = reshape(dmid,vmid.dim(1),vmid.dim(2),vmid.dim(3));
MATOUTTmix = zeros(vbg.dim);
for i = 1:size(indnewsub,1)
    indu = indnewsubT(i,:);
    induextendx = round(indu(1)-diffx):round(indu(1)+diffx);
    induextendy = round(indu(2)-diffy):round(indu(2)+diffy);
    induextendz = round(indu(3)-diffz):round(indu(3)+diffz);
    MATOUTTmix(induextendx,induextendy,induextendz) = dmidre(indnew(i));
end
dmidindex = unique(dmid);
dshownum = length(dmidindex)-1;
DATBG = dbgre;
DATBG(DATBG<30*0.9) = 0;
DATBG(DATBG>0&DATBG<=30*0.9) = 30+(dshownum)/(85-30);
DATBG(DATBG>=85) = 85-dshownum/55;
DATBG(DATBG>0) = (DATBG(DATBG>0)-30)/(85-30)*(dshownum);
DATBGwhole = DATBG;
DATBGwhole(MATOUTTmix>0) = MATOUTTmix(MATOUTTmix>0)-0.5+dshownum;
DATBGshow{1} = DATBGwhole;
listind{1,1} = 'Whole';

MAXMINX(1,:) = maxminx;
MAXMINY(1,:) = maxminy;
MAXMINZ(1,:) = maxminz;
if targetnum>1
%     MAXMINX(1,:) = maxminx;
%     MAXMINY(1,:) = maxminy;
%     MAXMINZ(1,:) = maxminz;
    for i = 1:targetnum
        listind{i+1,1} = ['TargetROI_',num2str(i)];
        dTARGET = zeros(size(MATOUTT));
        dTARGET(MATOUTT==dtarind(i+1)) = 1;
        MATOUTTmixtemp = MATOUTTmix.*dTARGET;
        MATOUTTmixsep{i} = MATOUTTmixtemp;
        DATBGwholetemp = DATBG;        
        DATBGwholetemp(MATOUTTmixtemp>0) = MATOUTTmixtemp(MATOUTTmixtemp>0)-0.5+dshownum;
        DATBGshow{i+1} = DATBGwholetemp;
        indt = find(dTARGET);
        [indtx indty indtz] = ind2sub(vbg.dim,indt);
        MAXMINX(i+1,:) = [min(indtx),max(indtx)];
        MAXMINY(i+1,:) = [min(indty),max(indty)];
        MAXMINZ(i+1,:) = [min(indtz),max(indtz)];
    end
else
    MATOUTTmixsep = [];
end

ASSHOWMAP.MAXMINX = MAXMINX;
ASSHOWMAP.MAXMINY = MAXMINY;
ASSHOWMAP.MAXMINZ = MAXMINZ;
save([ASSHOWMAP.Info.output,filesep,'Infoforshow',filesep,'Matrixshow.mat'],'DATBGshow','DATBG');
ASSHOWMAP.Showtarlist = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.05 0.4 0.9 0.45],'style','listbox','string',listind);
ASSHOWMAP.showind = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.05 0.35 0.9 0.04],'style','text','string',['sliceorder: min ',num2str(maxminz(1)),',max ',num2str(maxminz(2))]);
sliceorderdefault = maxminz(1):ceil((maxminz(2)-maxminz(1))/10):maxminz(2);
ASSHOWMAP.showinded = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.05 0.3 0.9 0.04],'style','edit','string',num2str(sliceorderdefault));
ASSHOWMAP.showbut = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.05 0.25 0.9 0.04],'style','pushbutton','string','show');
ASSHOWMAP.print = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.05 0.2 0.9 0.04],'style','pushbutton','string','print');
ASSHOWMAP.factortext = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.05 0.1 0.4 0.04],'style','text','string','Factor');
ASSHOWMAP.factoredit = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.55 0.1 0.4 0.04],'style','edit','string','0.2');
ASSHOWMAP.Return = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.05 0.05 0.9 0.04],'style','pushbutton','string','Return');
ASSHOWMAP.Matrix = uicontrol('parent',ASSHOWMAP.fig,'units','norm','pos',[0.05 0.15 0.9 0.04],'style','pushbutton','string','Matrix');
ASSHOWMAP.FigShow = figure('units','norm','pos',[0.12 0.07 0.8 0.8]);
ASSHOWMAP.FigShow2 = figure('units','norm','pos',[0.17 0.07 0.8 0.8]);
% sliceorderdefault
ASSHOWMAP.pat = pat;
% show default map
DATSHOWDEFAULT = DATBGshow{1};
load(fullfile(pat,'SEEDCOLOR.mat'));
load(fullfile(pat,'graycolmap.mat'));
indsep = round(1:63/dshownum:64);
colmaptemp = [];
for i = 1:dshownum
    colmaptemp(indsep(i):indsep(i+1),:) = ones(length(indsep(i):indsep(i+1)),1)*[SEEDCOLORSHOW(i,1),SEEDCOLORSHOW(i,2),SEEDCOLORSHOW(i,3)];
end
colmaptemp = [colgray;colmaptemp];
ASSHOWMAP.colmaptemp = colmaptemp;
asize = size(DATSHOWDEFAULT,2);
bsize = size(DATSHOWDEFAULT,1);
MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
for i = 1:length(sliceorderdefault)
    MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(DATSHOWDEFAULT(:,:,sliceorderdefault(i)));
    for j = 1:dshownum
        Showtempusedtemp = rot90(DATSHOWDEFAULT(:,:,sliceorderdefault(i)));
        Showtempusedtemp2 = rot90(DATBG(:,:,sliceorderdefault(i)));
        Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
        MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
    end
end
axestemp = axes('parent',ASSHOWMAP.FigShow,'units','norm','pos',[0.05 0.05 0.9 0.9]);
image(MATSHOWMATRIX,'parent',axestemp,'CDataMapping','scaled');axis(axestemp,'off');
colormap(axestemp,colmaptemp)
set(axestemp,'Clim',[0,dshownum*2])

POSX = MAXMINY(1,:);
POSY = MAXMINX(1,:);
LENX = abs(POSX(1)-POSX(2));
LENY = abs(POSY(1)-POSY(2));
FACTTEMP = get(ASSHOWMAP.factoredit,'string');
factor = str2num(FACTTEMP); % This is the change for the cut show.
FACTOR = 0.1;
ASSHOWMAP.factor = factor;
ASSHOWMAP.FACTOR = FACTOR;
sepfactorX = ceil(LENX*FACTOR);
sepfactorY = ceil(LENY*FACTOR);
halffactor = factor/2;
LENXh = ceil(LENX*halffactor);
LENYh = ceil(LENY*halffactor);
startpointX = POSX(1)-LENXh;
startpointY = POSY(1)-LENYh;
endpointX = POSX(2)+LENXh;
endpointY = POSY(2)+LENYh;
if startpointY<1
    startpointY = 1;
end
if startpointX<1
    startpointX = 1;
end
if endpointX>asize
    endpointX = asize;
end
if endpointY>bsize
    endpointY = bsize;
end
asizenew = length(startpointX:endpointX);
bsizenew = length(startpointY:endpointY);
MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
for i = 1:length(sliceorderdefault)
    MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
    for j = 1:dshownum
        Showtempusedtemp = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
        Showtempusedtemp2 = rot90(DATBG(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
        Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
        MATSHOWMATRIX2(1+asizenew*j+sepfactorX*j:asizenew*(j+1)+sepfactorX*j,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
    end
end
axestemp2 = axes('parent',ASSHOWMAP.FigShow2,'units','norm','pos',[0.05 0.05 0.9 0.9]);
image(MATSHOWMATRIX2,'parent',axestemp2,'CDataMapping','scaled');axis(axestemp2,'off');
colormap(axestemp2,colmaptemp)
set(axestemp2,'Clim',[0,dshownum*2])
%
ASSHOWMAP.dshownum = dshownum;
ASSHOWMAP.dmidindex = dmidindex;
%
%% callback
set(ASSHOWMAP.axi,'callback',{@AXIsel,ASSHOWMAP});
set(ASSHOWMAP.cor,'callback',{@CORsel,ASSHOWMAP});
set(ASSHOWMAP.sag,'callback',{@SAGsel,ASSHOWMAP});
set(ASSHOWMAP.showbut,'callback',{@Showmap,ASSHOWMAP});
set(ASSHOWMAP.print,'callback',{@Printmap,ASSHOWMAP});
set(ASSHOWMAP.Matrix,'callback',{@Matrixmap,ASSHOWMAP});
set(ASSHOWMAP.Showtarlist,'callback',{@Changesel,ASSHOWMAP});
set(ASSHOWMAP.Return,'callback',{@Returnmap,ASSHOWMAP});
end
function Radshow(varargin)
has = varargin{3};
INFODIR = get(has.infocoled,'string');
marks = find(INFODIR==':');
if ~isempty(marks)
    INFODIRS = INFODIR(marks+1:end);
else
    return;
end
load(INFODIRS);

[pat,nam,ext] = fileparts(which('AS_ShowWTAGUI.m'));
% bgpath = fullfile(pat,'mni_icbm152_t1_tal_nlin_asym_09a.nii');
% [vbg,dbg] = Dynamic_read_dir_NIFTI(bgpath);
% dbgre = reshape(dbg,vbg.dim(1),vbg.dim(2),vbg.dim(3));
% Info.bgpath = bgpath;
[vtar dtar] = Dynamic_read_dir_NIFTI(Info.target);
dtar(isnan(dtar)) = 0;
dtar(isinf(dtar)) = 0;
dtarind = unique(dtar);
dtarre = reshape(dtar,vtar.dim(1),vtar.dim(2),vtar.dim(3));
indnew = find(dtarre>0);
inmixdir = fullfile(Info.input,'mixedMaxID.nii');
[vmid dmid] = Dynamic_read_dir_NIFTI(inmixdir);
dmid(isnan(dmid)) = 0;
dmid(isinf(dmid)) = 0;
dmidindex = unique(dmid);
dshownum = length(dmidindex)-1;
for i = 1:dshownum
    Patternnum(i,1) = length(find(dmid==dmidindex(i+1)));
end
for i = 1:length(dtarind)-1
    for j = 1:dshownum
        PatternSub(i,j) = length(find(dmid(dtar==(dtarind(i+1)))==dmidindex(j+1)));
    end
end

Hsize = get(0,'screensize');
Patternperc = Patternnum./sum(Patternnum);
for i = 1:length(dtarind)-1
    PatternSubperc(i,:) = PatternSub(i,:)./sum(PatternSub(i,:));
end
load(fullfile(pat,'SEEDCOLOR.mat'));
msize = min(Hsize(3:4))*0.8;
H1 = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/4,msize,msize/2],'Name','Whole');
subplot(221);hold on;
for i = 1:dshownum
    bar(i,Patternnum(i),'barwidth',0.4,'facecolor',SEEDCOLORSHOW(i,:));
end

subplot(223);hold on;
for i = 1:dshownum
    bar(i,Patternperc(i),'barwidth',0.4,'facecolor',SEEDCOLORSHOW(i,:));
end
maxPATWHOLEPERC = max(Patternperc);
maxshow = ceil(maxPATWHOLEPERC*10)/10;
subplot(222);hold on;
for i = dshownum:-1:1
    barh(1,sum(Patternnum(1:i)),'barwidth',0.4,'facecolor',SEEDCOLORSHOW(i,:));
end
axis([0 sum(Patternnum)*1.1 0 2])
subplot(224);hold on;
for i = dshownum:-1:1
    barh(1,sum(Patternperc(1:i)),'barwidth',0.4,'facecolor',SEEDCOLORSHOW(i,:));
end
axis([0 1.1 0 2])

for i = 1:length(dtarind)-1
    HTemp{i} = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/4,msize,msize/2],'Name',['ROI: ', num2str(i)]);
    subplot(221);hold on;
    for j = 1:dshownum
        bar(j,PatternSub(i,j),'barwidth',0.4,'facecolor',SEEDCOLORSHOW(j,:));
    end
%     maxPATWHOLENUM = max(PatternSub);
    subplot(223);hold on;
    for j = 1:dshownum
        bar(j,PatternSubperc(i,j),'barwidth',0.4,'facecolor',SEEDCOLORSHOW(j,:));
    end
%     maxPATWHOLEPERC = max(Patternperc);
%     maxshow = ceil(maxPATWHOLEPERC*10)/10;
    subplot(222);hold on;
    for j = dshownum:-1:1
        barh(1,sum(PatternSub(i,1:j)),'barwidth',0.4,'facecolor',SEEDCOLORSHOW(j,:));
    end
    try
        axis([0 sum(PatternSub(i,:))*1.1 0 2])
    catch
        
    end
    subplot(224);hold on;
    for j = dshownum:-1:1
        barh(1,sum(PatternSubperc(i,1:j)),'barwidth',0.4,'facecolor',SEEDCOLORSHOW(j,:));
    end
    axis([0 1.1 0 2])
end

Hshowrad = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/2,msize,msize],'Name','Whole');
Hax = axes('parent',Hshowrad,'pos',[0.1 0.1 0.8 0.8]);
hold(Hax,'on')
maxPATWHOLENUM = max(Patternnum);
Outv = FindMaxNearVal(maxPATWHOLENUM);
% CirclePlot(Outv,Hax);
Radbg(Outv,Hax,dshownum)
PlotRadar(Patternnum,Hax,SEEDCOLORSHOW,dshownum)
axis(Hax,[-Outv,Outv,-Outv,Outv]*1.1)
Hshowrad2 = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/2,msize,msize],'Name','Whole');
Hax2 = axes('parent',Hshowrad2,'pos',[0.1 0.1 0.8 0.8]);
hold(Hax2,'on')
maxPATWHOLENUM = max(Patternnum);
Outv = FindMaxNearVal(maxPATWHOLENUM);
% CirclePlot(Outv,Hax);
Radbg2(Outv,Hax2,dshownum)
PlotRadar(Patternnum,Hax2,SEEDCOLORSHOW,dshownum)
axis(Hax2,[-Outv,Outv,-Outv,Outv]*1.1)
for i = 1:length(dtarind)-1
    Hshowradsep(i) = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/2,msize,msize],'Name',['ROI: ',num2str(i)]);
    Hax = axes('parent',Hshowradsep(i),'pos',[0.1 0.1 0.8 0.8]);
    hold(Hax,'on')
    TEMPV = PatternSub(i,:);
    maxPATWHOLENUM = max(TEMPV);
    Outv = FindMaxNearVal(maxPATWHOLENUM);
    % CirclePlot(Outv,Hax);
    Radbg(Outv,Hax,dshownum)
    PlotRadar(TEMPV,Hax,SEEDCOLORSHOW,dshownum)
    axis(Hax,[-Outv,Outv,-Outv,Outv]*1.1)
    Hshowradsep2(i) = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/2,msize,msize],'Name',['ROI: ',num2str(i)]);
    Hax2 = axes('parent',Hshowradsep2(i),'pos',[0.1 0.1 0.8 0.8]);
    hold(Hax2,'on')
    maxPATWHOLENUM = max(TEMPV);
    Outv = FindMaxNearVal(maxPATWHOLENUM);
    % CirclePlot(Outv,Hax);
    Radbg2(Outv,Hax2,dshownum)
    PlotRadar(TEMPV,Hax2,SEEDCOLORSHOW,dshownum)
    axis(Hax2,[-Outv,Outv,-Outv,Outv]*1.1)
end
%%

Hshowrad = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/2,msize,msize],'Name','Percent Whole');
Hax = axes('parent',Hshowrad,'pos',[0.1 0.1 0.8 0.8]);
hold(Hax,'on')
maxPATWHOLENUM = max(Patternperc);
Outv = FindMaxNearVal(maxPATWHOLENUM);
% Outv = FindMaxNearValper(maxPATWHOLENUM);
% CirclePlot(Outv,Hax);
Radbg(Outv,Hax,dshownum)
PlotRadar(Patternperc,Hax,SEEDCOLORSHOW,dshownum)
axis(Hax,[-Outv,Outv,-Outv,Outv]*1.1)
Hshowrad2 = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/2,msize,msize],'Name','Percent Whole');
Hax2 = axes('parent',Hshowrad2,'pos',[0.1 0.1 0.8 0.8]);
hold(Hax2,'on')
maxPATWHOLENUM = max(Patternperc);
Outv = FindMaxNearVal(maxPATWHOLENUM);
% Outv = FindMaxNearValper(maxPATWHOLENUM);
% CirclePlot(Outv,Hax);
Radbg2(Outv,Hax2,dshownum)
PlotRadar(Patternperc,Hax2,SEEDCOLORSHOW,dshownum)
axis(Hax2,[-Outv,Outv,-Outv,Outv]*1.1)
for i = 1:length(dtarind)-1
    Hshowradsep(i) = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/2,msize,msize],'Name',['Percent ROI: ',num2str(i)]);
    Hax = axes('parent',Hshowradsep(i),'pos',[0.1 0.1 0.8 0.8]);
    hold(Hax,'on')
    TEMPV = PatternSubperc(i,:);
    maxPATWHOLENUM = max(TEMPV);
    Outv = FindMaxNearVal(maxPATWHOLENUM);
%     Outv = FindMaxNearValper(maxPATWHOLENUM);
    % CirclePlot(Outv,Hax);
    Radbg(Outv,Hax,dshownum)
    PlotRadar(TEMPV,Hax,SEEDCOLORSHOW,dshownum)
    axis(Hax,[-Outv,Outv,-Outv,Outv]*1.1)
    Hshowradsep2(i) = figure('pos',[Hsize(3)/2-msize/2,Hsize(4)/2-msize/2,msize,msize],'Name',['Percent ROI: ',num2str(i)]);
    Hax2 = axes('parent',Hshowradsep2(i),'pos',[0.1 0.1 0.8 0.8]);
    hold(Hax2,'on')
    maxPATWHOLENUM = max(TEMPV);
    Outv = FindMaxNearVal(maxPATWHOLENUM);
%     Outv = FindMaxNearValper(maxPATWHOLENUM);
    % CirclePlot(Outv,Hax);
    Radbg2(Outv,Hax2,dshownum)
    PlotRadar(TEMPV,Hax2,SEEDCOLORSHOW,dshownum)
    axis(Hax2,[-Outv,Outv,-Outv,Outv]*1.1)
end
end
function Exithas(varargin)
has = varargin{3};
close(has.fig)
BCCT_VIEW;
end
%%
function Returnmap(varargin)
ASSHOWMAP = varargin{3};
close(ASSHOWMAP.fig);
close(ASSHOWMAP.FigShow)
close(ASSHOWMAP.FigShow2);
end
function Changesel(varargin)
ASSHOWMAP = varargin{3};
targetnum = ASSHOWMAP.Info.targetnum;
if targetnum>1
    Valuse = get(ASSHOWMAP.Showtarlist,'val');
    Valaxi = get(ASSHOWMAP.axi,'val');
    Valcor = get(ASSHOWMAP.cor,'val');
    Valsag = get(ASSHOWMAP.sag,'val');
    if Valaxi
        set(ASSHOWMAP.showind,'string',['sliceorder: min ',num2str(ASSHOWMAP.MAXMINZ(Valuse,1)),',max ',num2str(ASSHOWMAP.MAXMINZ(Valuse,2))])
        MMZ = ASSHOWMAP.MAXMINZ(Valuse,:);
        sliceorderdefault = MMZ(1):ceil((MMZ(2)-MMZ(1))/10):MMZ(2);
        set(ASSHOWMAP.showinded,'string',num2str(sliceorderdefault));
    elseif Valcor
        set(ASSHOWMAP.showind,'string',['sliceorder: min ',num2str(ASSHOWMAP.MAXMINY(Valuse,1)),',max ',num2str(ASSHOWMAP.MAXMINY(Valuse,2))])
        MMY = ASSHOWMAP.MAXMINY(Valuse,:);  
        sliceorderdefault = MMY(1):ceil((MMY(2)-MMY(1))/10):MMY(2);
        set(ASSHOWMAP.showinded,'string',num2str(sliceorderdefault));
    elseif Valsag
        set(ASSHOWMAP.showind,'string',['sliceorder: min ',num2str(ASSHOWMAP.MAXMINX(Valuse,1)),',max ',num2str(ASSHOWMAP.MAXMINX(Valuse,2))]) 
        MMX = ASSHOWMAP.MAXMINX(Valuse,:);
        sliceorderdefault = MMX(1):ceil((MMX(2)-MMX(1))/10):MMX(2);    
        set(ASSHOWMAP.showinded,'string',num2str(sliceorderdefault));     
    end
end
end
function AXIsel(varargin)
ASSHOWMAP = varargin{3};
set(ASSHOWMAP.axi,'val',1);
set(ASSHOWMAP.cor,'val',0);
set(ASSHOWMAP.sag,'val',0);
MAXMINX = ASSHOWMAP.MAXMINX;
MAXMINY = ASSHOWMAP.MAXMINY;
MAXMINZ = ASSHOWMAP.MAXMINZ;
colmaptemp = ASSHOWMAP.colmaptemp;
load([ASSHOWMAP.Info.output,filesep,'Infoforshow',filesep,'Matrixshow.mat']);
dshownum = ASSHOWMAP.dshownum;
dmidindex = ASSHOWMAP.dmidindex;
sliceorderdefault = MAXMINZ(1,1):ceil((MAXMINZ(1,2)-MAXMINZ(1,1))/10):MAXMINZ(1,2);
set(ASSHOWMAP.showind,'string',['sliceorder: min ',num2str(MAXMINZ(1,1)),',max ',num2str(MAXMINZ(1,2))]);
set(ASSHOWMAP.showinded,'string',num2str(sliceorderdefault));
set(ASSHOWMAP.Showtarlist,'val',1);
DATSHOWDEFAULT = DATBGshow{1};
asize = size(DATSHOWDEFAULT,2);
bsize = size(DATSHOWDEFAULT,1);
MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
for i = 1:length(sliceorderdefault)
    MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(DATSHOWDEFAULT(:,:,sliceorderdefault(i)));
    for j = 1:dshownum
        Showtempusedtemp = rot90(DATSHOWDEFAULT(:,:,sliceorderdefault(i)));
        Showtempusedtemp2 = rot90(DATBG(:,:,sliceorderdefault(i)));
        Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
        MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
    end
end
axestemp = axes('parent',ASSHOWMAP.FigShow,'units','norm','pos',[0.05 0.05 0.9 0.9]);
image(MATSHOWMATRIX,'parent',axestemp,'CDataMapping','scaled');axis(axestemp,'off');
colormap(axestemp,colmaptemp)
set(axestemp,'Clim',[0,dshownum*2])

POSX = MAXMINY(1,:);
POSY = MAXMINX(1,:);
LENX = abs(POSX(1)-POSX(2));
LENY = abs(POSY(1)-POSY(2));
FACTTEMP = get(ASSHOWMAP.factoredit,'string');
factor = str2num(FACTTEMP); % This is the change for the cut show.
FACTOR = 0.1;
%     ASSHOWMAP.factor = factor;
%     ASSHOWMAP.FACTOR = FACTOR;
sepfactorX = ceil(LENX*FACTOR);
sepfactorY = ceil(LENY*FACTOR);
halffactor = factor/2;
LENXh = ceil(LENX*halffactor);
LENYh = ceil(LENY*halffactor);
startpointX = POSX(1)-LENXh;
startpointY = POSY(1)-LENYh;
endpointX = POSX(2)+LENXh;
endpointY = POSY(2)+LENYh;
if startpointY<1
    startpointY = 1;
end
if startpointX<1
    startpointX = 1;
end
if endpointX>asize
    endpointX = asize;
end
if endpointY>bsize
    endpointY = bsize;
end
asizenew = length(startpointX:endpointX);
bsizenew = length(startpointY:endpointY);
MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
for i = 1:length(sliceorderdefault)
    MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
    for j = 1:dshownum
        Showtempusedtemp = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
        Showtempusedtemp2 = rot90(DATBG(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
        Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
        MATSHOWMATRIX2(1+asizenew*j+sepfactorX*(j+1):asizenew*(j+1)+sepfactorX*(j+1),1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
    end
end
axestemp2 = axes('parent',ASSHOWMAP.FigShow2,'units','norm','pos',[0.05 0.05 0.9 0.9]);
image(MATSHOWMATRIX2,'parent',axestemp2,'CDataMapping','scaled');axis(axestemp2,'off');
colormap(axestemp2,colmaptemp)
set(axestemp2,'Clim',[0,dshownum*2])
end
function CORsel(varargin)
ASSHOWMAP = varargin{3};
set(ASSHOWMAP.axi,'val',0);
set(ASSHOWMAP.cor,'val',1);
set(ASSHOWMAP.sag,'val',0);
MAXMINX = ASSHOWMAP.MAXMINX;
MAXMINY = ASSHOWMAP.MAXMINY;
MAXMINZ = ASSHOWMAP.MAXMINZ;
colmaptemp = ASSHOWMAP.colmaptemp;
load([ASSHOWMAP.Info.output,filesep,'Infoforshow',filesep,'Matrixshow.mat']);
dshownum = ASSHOWMAP.dshownum;
dmidindex = ASSHOWMAP.dmidindex;
sliceorderdefault = MAXMINY(1,1):ceil((MAXMINY(1,2)-MAXMINY(1,1))/10):MAXMINY(1,2);
set(ASSHOWMAP.showind,'string',['sliceorder: min ',num2str(MAXMINY(1,1)),',max ',num2str(MAXMINY(1,2))]);
set(ASSHOWMAP.showinded,'string',num2str(sliceorderdefault));
set(ASSHOWMAP.Showtarlist,'val',1);
DATSHOWDEFAULT = DATBGshow{1};
asize = size(DATSHOWDEFAULT,3);
bsize = size(DATSHOWDEFAULT,1);
MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
for i = 1:length(sliceorderdefault)
    MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(squeeze(DATSHOWDEFAULT(:,sliceorderdefault(i),:)));
    for j = 1:dshownum
        Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(:,sliceorderdefault(i),:)));
        Showtempusedtemp2 = rot90(squeeze(DATBG(:,sliceorderdefault(i),:)));
        Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
        MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
    end
end
axestemp = axes('parent',ASSHOWMAP.FigShow,'units','norm','pos',[0.05 0.05 0.9 0.9]);
image(MATSHOWMATRIX,'parent',axestemp,'CDataMapping','scaled');axis(axestemp,'off');
colormap(axestemp,colmaptemp)
set(axestemp,'Clim',[0,dshownum*2])

POSX = MAXMINZ(1,:);
POSY = MAXMINX(1,:);
LENX = abs(POSX(1)-POSX(2));
LENY = abs(POSY(1)-POSY(2));
FACTTEMP = get(ASSHOWMAP.factoredit,'string');
factor = str2num(FACTTEMP); % This is the change for the cut show.
FACTOR = 0.1;
%     ASSHOWMAP.factor = factor;
%     ASSHOWMAP.FACTOR = FACTOR;
sepfactorX = ceil(LENX*FACTOR);
sepfactorY = ceil(LENY*FACTOR);
halffactor = factor/2;
LENXh = ceil(LENX*halffactor);
LENYh = ceil(LENY*halffactor);
startpointX = POSX(1)-LENXh;
startpointY = POSY(1)-LENYh;
endpointX = POSX(2)+LENXh;
endpointY = POSY(2)+LENYh;
if startpointY<1
    startpointY = 1;
end
if startpointX<1
    startpointX = 1;
end
if endpointX>asize
    endpointX = asize;
end
if endpointY>bsize
    endpointY = bsize;
end
asizenew = length(startpointX:endpointX);
bsizenew = length(startpointY:endpointY);
MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
for i = 1:length(sliceorderdefault)
    MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(squeeze(DATSHOWDEFAULT(startpointY:endpointY,sliceorderdefault(i),startpointX:endpointX)));
    for j = 1:dshownum
        Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(startpointY:endpointY,sliceorderdefault(i),startpointX:endpointX)));
        Showtempusedtemp2 = rot90(squeeze(DATBG(startpointY:endpointY,sliceorderdefault(i),startpointX:endpointX)));
        Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
        MATSHOWMATRIX2(1+asizenew*j+sepfactorX*(j+1):asizenew*(j+1)+sepfactorX*(j+1),1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
    end
end
axestemp2 = axes('parent',ASSHOWMAP.FigShow2,'units','norm','pos',[0.05 0.05 0.9 0.9]);
image(MATSHOWMATRIX2,'parent',axestemp2,'CDataMapping','scaled');axis(axestemp2,'off');
colormap(axestemp2,colmaptemp)
set(axestemp2,'Clim',[0,dshownum*2])

end
function SAGsel(varargin)
ASSHOWMAP = varargin{3};
set(ASSHOWMAP.axi,'val',0);
set(ASSHOWMAP.cor,'val',0);
set(ASSHOWMAP.sag,'val',1);
MAXMINX = ASSHOWMAP.MAXMINX;
MAXMINY = ASSHOWMAP.MAXMINY;
MAXMINZ = ASSHOWMAP.MAXMINZ;
colmaptemp = ASSHOWMAP.colmaptemp;
load([ASSHOWMAP.Info.output,filesep,'Infoforshow',filesep,'Matrixshow.mat']);
dshownum = ASSHOWMAP.dshownum;
dmidindex = ASSHOWMAP.dmidindex;
sliceorderdefault = MAXMINX(1,1):ceil((MAXMINX(1,2)-MAXMINX(1,1))/10):MAXMINX(1,2);
set(ASSHOWMAP.showind,'string',['sliceorder: min ',num2str(MAXMINX(1,1)),',max ',num2str(MAXMINX(1,2))]);
set(ASSHOWMAP.showinded,'string',num2str(sliceorderdefault));
set(ASSHOWMAP.Showtarlist,'val',1);
DATSHOWDEFAULT = DATBGshow{1};
asize = size(DATSHOWDEFAULT,3);
bsize = size(DATSHOWDEFAULT,2);
MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
for i = 1:length(sliceorderdefault)
    MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),:,:)));
    for j = 1:dshownum
        Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),:,:)));
        Showtempusedtemp2 = rot90(squeeze(DATBG(sliceorderdefault(i),:,:)));
        Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
        MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
    end
end
axestemp = axes('parent',ASSHOWMAP.FigShow,'units','norm','pos',[0.05 0.05 0.9 0.9]);
image(MATSHOWMATRIX,'parent',axestemp,'CDataMapping','scaled');axis(axestemp,'off');
colormap(axestemp,colmaptemp)
set(axestemp,'Clim',[0,dshownum*2])

POSX = MAXMINZ(1,:);
POSY = MAXMINY(1,:);
LENX = abs(POSX(1)-POSX(2));
LENY = abs(POSY(1)-POSY(2));
FACTTEMP = get(ASSHOWMAP.factoredit,'string');
factor = str2num(FACTTEMP); % This is the change for the cut show.
FACTOR = 0.1;
ASSHOWMAP.factor = factor;
ASSHOWMAP.FACTOR = FACTOR;
sepfactorX = ceil(LENX*FACTOR);
sepfactorY = ceil(LENY*FACTOR);
halffactor = factor/2;
LENXh = ceil(LENX*halffactor);
LENYh = ceil(LENY*halffactor);
startpointX = POSX(1)-LENXh;
startpointY = POSY(1)-LENYh;
endpointX = POSX(2)+LENXh;
endpointY = POSY(2)+LENYh;
if startpointY<1
    startpointY = 1;
end
if startpointX<1
    startpointX = 1;
end
if endpointX>asize
    endpointX = asize;
end
if endpointY>bsize
    endpointY = bsize;
end
asizenew = length(startpointX:endpointX);
bsizenew = length(startpointY:endpointY);
MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
for i = 1:length(sliceorderdefault)
    MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),startpointY:endpointY,startpointX:endpointX)));
    for j = 1:dshownum
        Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),startpointY:endpointY,startpointX:endpointX)));
        Showtempusedtemp2 = rot90(squeeze(DATBG(sliceorderdefault(i),startpointY:endpointY,startpointX:endpointX)));
        Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
        MATSHOWMATRIX2(1+asizenew*j+sepfactorX*(j+1):asizenew*(j+1)+sepfactorX*(j+1),1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
    end
end
axestemp2 = axes('parent',ASSHOWMAP.FigShow2,'units','norm','pos',[0.05 0.05 0.9 0.9]);
image(MATSHOWMATRIX2,'parent',axestemp2,'CDataMapping','scaled');axis(axestemp2,'off');
colormap(axestemp2,colmaptemp)
set(axestemp2,'Clim',[0,dshownum*2])
end
function Showmap(varargin)
ASSHOWMAP = varargin{3};
axiv = get(ASSHOWMAP.axi,'val');
corv = get(ASSHOWMAP.cor,'val');
sagv = get(ASSHOWMAP.sag,'val');
lists = get(ASSHOWMAP.Showtarlist,'val');
targetnum = ASSHOWMAP.Info.targetnum;

sliceinfomations = get(ASSHOWMAP.showinded,'string');
sliceorderdefault = str2num(sliceinfomations);

% save([ASSHOWMAP.Info.output,filesep,'Infoforshow',filesep,'Matrixshow.mat'],'DATBGshow','DATBG');
load([ASSHOWMAP.Info.output,filesep,'Infoforshow',filesep,'Matrixshow.mat']);
dshownum = ASSHOWMAP.dshownum;
dmidindex = ASSHOWMAP.dmidindex;
MAXMINX = ASSHOWMAP.MAXMINX;
MAXMINY = ASSHOWMAP.MAXMINY;
MAXMINZ = ASSHOWMAP.MAXMINZ;
colmaptemp = ASSHOWMAP.colmaptemp;
%%

DATSHOWDEFAULT = DATBGshow{lists};
if axiv
    asize = size(DATSHOWDEFAULT,2);
    bsize = size(DATSHOWDEFAULT,1);
    MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(DATSHOWDEFAULT(:,:,sliceorderdefault(i)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(DATSHOWDEFAULT(:,:,sliceorderdefault(i)));
            Showtempusedtemp2 = rot90(DATBG(:,:,sliceorderdefault(i)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
        end
    end
    axestemp = axes('parent',ASSHOWMAP.FigShow,'units','norm','pos',[0.05 0.05 0.9 0.9]);
    image(MATSHOWMATRIX,'parent',axestemp,'CDataMapping','scaled');axis(axestemp,'off');
    colormap(axestemp,colmaptemp)
    set(axestemp,'Clim',[0,dshownum*2])
    
    POSX = MAXMINY(lists,:);
    POSY = MAXMINX(lists,:);
    LENX = abs(POSX(1)-POSX(2));
    LENY = abs(POSY(1)-POSY(2));
    FACTTEMP = get(ASSHOWMAP.factoredit,'string');
    factor = str2num(FACTTEMP); % This is the change for the cut show.
    FACTOR = 0.1;
%     ASSHOWMAP.factor = factor;
%     ASSHOWMAP.FACTOR = FACTOR;
    sepfactorX = ceil(LENX*FACTOR);
    sepfactorY = ceil(LENY*FACTOR);
    halffactor = factor/2;
    LENXh = ceil(LENX*halffactor);
    LENYh = ceil(LENY*halffactor);
    startpointX = POSX(1)-LENXh;
    startpointY = POSY(1)-LENYh;
    endpointX = POSX(2)+LENXh;
    endpointY = POSY(2)+LENYh;
    if startpointY<1
        startpointY = 1;
    end
    if startpointX<1
        startpointX = 1;
    end
    if endpointX>asize
        endpointX = asize;
    end
    if endpointY>bsize
        endpointY = bsize;
    end
    asizenew = length(startpointX:endpointX);
    bsizenew = length(startpointY:endpointY);
    MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
            Showtempusedtemp2 = rot90(DATBG(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX2(1+asizenew*j+sepfactorX*(j+1):asizenew*(j+1)+sepfactorX*(j+1),1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
        end
    end
    axestemp2 = axes('parent',ASSHOWMAP.FigShow2,'units','norm','pos',[0.05 0.05 0.9 0.9]);
    image(MATSHOWMATRIX2,'parent',axestemp2,'CDataMapping','scaled');axis(axestemp2,'off');
    colormap(axestemp2,colmaptemp)
    set(axestemp2,'Clim',[0,dshownum*2])
elseif corv
    asize = size(DATSHOWDEFAULT,3);
    bsize = size(DATSHOWDEFAULT,1);
    MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(squeeze(DATSHOWDEFAULT(:,sliceorderdefault(i),:)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(:,sliceorderdefault(i),:)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(:,sliceorderdefault(i),:)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
        end
    end
    axestemp = axes('parent',ASSHOWMAP.FigShow,'units','norm','pos',[0.05 0.05 0.9 0.9]);
    image(MATSHOWMATRIX,'parent',axestemp,'CDataMapping','scaled');axis(axestemp,'off');
    colormap(axestemp,colmaptemp)
    set(axestemp,'Clim',[0,dshownum*2])
    
    POSX = MAXMINZ(lists,:);
    POSY = MAXMINX(lists,:);
    LENX = abs(POSX(1)-POSX(2));
    LENY = abs(POSY(1)-POSY(2));
    FACTTEMP = get(ASSHOWMAP.factoredit,'string');
    factor = str2num(FACTTEMP); % This is the change for the cut show.
    FACTOR = 0.1;
%     ASSHOWMAP.factor = factor;
%     ASSHOWMAP.FACTOR = FACTOR;
    sepfactorX = ceil(LENX*FACTOR);
    sepfactorY = ceil(LENY*FACTOR);
    halffactor = factor/2;
    LENXh = ceil(LENX*halffactor);
    LENYh = ceil(LENY*halffactor);
    startpointX = POSX(1)-LENXh;
    startpointY = POSY(1)-LENYh;
    endpointX = POSX(2)+LENXh;
    endpointY = POSY(2)+LENYh;
    if startpointY<1
        startpointY = 1;
    end
    if startpointX<1
        startpointX = 1;
    end
    if endpointX>asize
        endpointX = asize;
    end
    if endpointY>bsize
        endpointY = bsize;
    end
    asizenew = length(startpointX:endpointX);
    bsizenew = length(startpointY:endpointY);
    MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(squeeze(DATSHOWDEFAULT(startpointY:endpointY,sliceorderdefault(i),startpointX:endpointX)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(startpointY:endpointY,sliceorderdefault(i),startpointX:endpointX)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(startpointY:endpointY,sliceorderdefault(i),startpointX:endpointX)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX2(1+asizenew*j+sepfactorX*(j+1):asizenew*(j+1)+sepfactorX*(j+1),1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
        end
    end
    axestemp2 = axes('parent',ASSHOWMAP.FigShow2,'units','norm','pos',[0.05 0.05 0.9 0.9]);
    image(MATSHOWMATRIX2,'parent',axestemp2,'CDataMapping','scaled');axis(axestemp2,'off');
    colormap(axestemp2,colmaptemp)
    set(axestemp2,'Clim',[0,dshownum*2])    
elseif sagv
    asize = size(DATSHOWDEFAULT,3);
    bsize = size(DATSHOWDEFAULT,2);
    MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),:,:)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),:,:)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(sliceorderdefault(i),:,:)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
        end
    end
    axestemp = axes('parent',ASSHOWMAP.FigShow,'units','norm','pos',[0.05 0.05 0.9 0.9]);
    image(MATSHOWMATRIX,'parent',axestemp,'CDataMapping','scaled');axis(axestemp,'off');
    colormap(axestemp,colmaptemp)
    set(axestemp,'Clim',[0,dshownum*2])
    
    POSX = MAXMINZ(lists,:);
    POSY = MAXMINY(lists,:);
    LENX = abs(POSX(1)-POSX(2));
    LENY = abs(POSY(1)-POSY(2));
    FACTTEMP = get(ASSHOWMAP.factoredit,'string');
    factor = str2num(FACTTEMP); % This is the change for the cut show.
    FACTOR = 0.1;
    ASSHOWMAP.factor = factor;
    ASSHOWMAP.FACTOR = FACTOR;
    sepfactorX = ceil(LENX*FACTOR);
    sepfactorY = ceil(LENY*FACTOR);
    halffactor = factor/2;
    LENXh = ceil(LENX*halffactor);
    LENYh = ceil(LENY*halffactor);
    startpointX = POSX(1)-LENXh;
    startpointY = POSY(1)-LENYh;
    endpointX = POSX(2)+LENXh;
    endpointY = POSY(2)+LENYh;
    if startpointY<1
        startpointY = 1;
    end
    if startpointX<1
        startpointX = 1;
    end
    if endpointX>asize
        endpointX = asize;
    end
    if endpointY>bsize
        endpointY = bsize;
    end
    asizenew = length(startpointX:endpointX);
    bsizenew = length(startpointY:endpointY);
    MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),startpointY:endpointY,startpointX:endpointX)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),startpointY:endpointY,startpointX:endpointX)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(sliceorderdefault(i),startpointY:endpointY,startpointX:endpointX)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX2(1+asizenew*j+sepfactorX*(j+1):asizenew*(j+1)+sepfactorX*(j+1),1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
        end
    end
    axestemp2 = axes('parent',ASSHOWMAP.FigShow2,'units','norm','pos',[0.05 0.05 0.9 0.9]);
    image(MATSHOWMATRIX2,'parent',axestemp2,'CDataMapping','scaled');axis(axestemp2,'off');
    colormap(axestemp2,colmaptemp)
    set(axestemp2,'Clim',[0,dshownum*2])    
end

end
function Printmap(varargin)
ASSHOWMAP = varargin{3};
axiv = get(ASSHOWMAP.axi,'val');
corv = get(ASSHOWMAP.cor,'val');
sagv = get(ASSHOWMAP.sag,'val');
lists = get(ASSHOWMAP.Showtarlist,'val');
targetnum = ASSHOWMAP.Info.targetnum;

sliceinfomations = get(ASSHOWMAP.showinded,'string');
sliceorderdefault = str2num(sliceinfomations);
% save testinfo
% save([ASSHOWMAP.Info.output,filesep,'Infoforshow',filesep,'Matrixshow.mat'],'DATBGshow','DATBG');
load([ASSHOWMAP.Info.output,filesep,'Infoforshow',filesep,'Matrixshow.mat']);
dshownum = ASSHOWMAP.dshownum;
dmidindex = ASSHOWMAP.dmidindex;
MAXMINX = ASSHOWMAP.MAXMINX;
MAXMINY = ASSHOWMAP.MAXMINY;
MAXMINZ = ASSHOWMAP.MAXMINZ;
colmaptemp = ASSHOWMAP.colmaptemp;
OUTDIRS = ASSHOWMAP.Info.output;
OUTDIRSDIR1 = [OUTDIRS,filesep,'Axial'];
OUTDIRSDIR2 = [OUTDIRS,filesep,'Cornoral'];
OUTDIRSDIR3 = [OUTDIRS,filesep,'Sagittal'];
mkdir(OUTDIRSDIR1);
mkdir(OUTDIRSDIR2);
mkdir(OUTDIRSDIR3);
Hsize = get(0,'screensize');
Hexist1 = Hsize(3)-200;
Hexist2 = Hsize(4)-200;
%%

DATSHOWDEFAULT = DATBGshow{lists};
if axiv
    asize = size(DATSHOWDEFAULT,2);
    bsize = size(DATSHOWDEFAULT,1);
    MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(DATSHOWDEFAULT(:,:,sliceorderdefault(i)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(DATSHOWDEFAULT(:,:,sliceorderdefault(i)));
            Showtempusedtemp2 = rot90(DATBG(:,:,sliceorderdefault(i)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
        end
    end
    Dsize = [size(MATSHOWMATRIX,2),size(MATSHOWMATRIX,1)];
    factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
    H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
    imagesc(MATSHOWMATRIX,[0,dshownum*2]);colormap(colmaptemp);
    axis off;
    saveas(H,[OUTDIRSDIR1,filesep,'CombinedWhole.fig'])
    set(H,'PaperPositionMode','manual');
    set(H,'PaperUnits','inch')
    XSIZE = size(MATSHOWMATRIX,2);
    YSIZE = size(MATSHOWMATRIX,1);
    factor = 1:100;
    XSIZEnew = XSIZE*factor;
    YSIZEnew = YSIZE*factor;
    FACTORS1 = find(XSIZEnew>Hsize(3));
    FACTORS2 = find(YSIZEnew>Hsize(4));
    FACTORS = max(FACTORS1(1),FACTORS2(1));
    XSIZEU = XSIZEnew(FACTORS);
    YSIZEU = YSIZEnew(FACTORS);
    set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
    print(H,[OUTDIRSDIR1,filesep,'CombinedWhole.tif'],'-dtiff','-r300')
    close(H)
    %%
    SingleIND = MAXMINZ(lists,1):MAXMINZ(lists,2);
    for i = SingleIND
        asize = size(DATSHOWDEFAULT,2);
        bsize = size(DATSHOWDEFAULT,1);        
        TempOut = rot90(DATSHOWDEFAULT(:,:,i));
        Dsize = [size(TempOut,2),size(TempOut,1)];
        factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
        H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
        imagesc(TempOut,[0,dshownum*2]);colormap(colmaptemp);
        axis off;
        saveas(H,[OUTDIRSDIR1,filesep,'SepWholeSlice',num2str(i),'.fig'])
        set(H,'PaperPositionMode','manual');
        set(H,'PaperUnits','inch')
        XSIZE = size(TempOut,2);
        YSIZE = size(TempOut,1);
        factor = 1:100;
        XSIZEnew = XSIZE*factor;
        YSIZEnew = YSIZE*factor;
        FACTORS1 = find(XSIZEnew>Hsize(3));
        FACTORS2 = find(YSIZEnew>Hsize(4));
        FACTORS = max(FACTORS1(1),FACTORS2(1));
        XSIZEU = XSIZEnew(FACTORS);
        YSIZEU = YSIZEnew(FACTORS);
        set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
        print(H,[OUTDIRSDIR1,filesep,'SepWholeSlice',num2str(i),'.tif'],'-dtiff','-r300')
        close(H)
        for j = 1:dshownum
            Showtempusedtemp = rot90(DATSHOWDEFAULT(:,:,i));
            Showtempusedtemp2 = rot90(DATBG(:,:,i));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            Dsize = [size(Showtempusedtemp2,2),size(Showtempusedtemp2,1)];
            factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
            H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
            imagesc(Showtempusedtemp2,[0,dshownum*2]);colormap(colmaptemp);
            axis off;
            saveas(H,[OUTDIRSDIR1,filesep,'SepWholeSlice',num2str(i),'ROI',num2str(j),'.fig'])
            set(H,'PaperPositionMode','manual');
            set(H,'PaperUnits','inch')
            XSIZE = size(Showtempusedtemp2,2);
            YSIZE = size(Showtempusedtemp2,1);
            factor = 1:100;
            XSIZEnew = XSIZE*factor;
            YSIZEnew = YSIZE*factor;
            FACTORS1 = find(XSIZEnew>Hsize(3));
            FACTORS2 = find(YSIZEnew>Hsize(4));
            FACTORS = max(FACTORS1(1),FACTORS2(1));
            XSIZEU = XSIZEnew(FACTORS);
            YSIZEU = YSIZEnew(FACTORS);
            set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
            print(H,[OUTDIRSDIR1,filesep,'SepWholeSlice',num2str(i),'ROI',num2str(j),'.tif'],'-dtiff','-r300')
            close(H)
        end
    end
    %%
    POSX = MAXMINY(lists,:);
    POSY = MAXMINX(lists,:);
    LENX = abs(POSX(1)-POSX(2));
    LENY = abs(POSY(1)-POSY(2));
    FACTTEMP = get(ASSHOWMAP.factoredit,'string');
    factor = str2num(FACTTEMP); % This is the change for the cut show.
    FACTOR = 0.1;
    sepfactorX = ceil(LENX*FACTOR);
    sepfactorY = ceil(LENY*FACTOR);
    halffactor = factor/2;
    LENXh = ceil(LENX*halffactor);
    LENYh = ceil(LENY*halffactor);
    startpointX = POSX(1)-LENXh;
    startpointY = POSY(1)-LENYh;
    endpointX = POSX(2)+LENXh;
    endpointY = POSY(2)+LENYh;
    if startpointY<1
        startpointY = 1;
    end
    if startpointX<1
        startpointX = 1;
    end
    if endpointX>asize
        endpointX = asize;
    end
    if endpointY>bsize
        endpointY = bsize;
    end
    asizenew = length(startpointX:endpointX);
    bsizenew = length(startpointY:endpointY);
    MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
            Showtempusedtemp2 = rot90(DATBG(startpointY:endpointY,startpointX:endpointX,sliceorderdefault(i)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX2(1+asizenew*j+sepfactorX*(j+1):asizenew*(j+1)+sepfactorX*(j+1),1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
        end
    end
    Dsize = [size(MATSHOWMATRIX2,2),size(MATSHOWMATRIX2,1)];
    factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
    H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
    imagesc(MATSHOWMATRIX2,[0,dshownum*2]);colormap(colmaptemp);
    axis off;
    saveas(H,[OUTDIRSDIR1,filesep,'CombinedRegion.fig'])
    set(H,'PaperPositionMode','manual');
    set(H,'PaperUnits','inch')
    XSIZE = size(MATSHOWMATRIX2,2);
    YSIZE = size(MATSHOWMATRIX2,1);
    factor = 1:100;
    XSIZEnew = XSIZE*factor;
    YSIZEnew = YSIZE*factor;
    FACTORS1 = find(XSIZEnew>Hsize(3));
    FACTORS2 = find(YSIZEnew>Hsize(4));
    FACTORS = max(FACTORS1(1),FACTORS2(1));
    XSIZEU = XSIZEnew(FACTORS);
    YSIZEU = YSIZEnew(FACTORS);
    set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
    print(H,[OUTDIRSDIR1,filesep,'CombinedRegion.tif'],'-dtiff','-r300')
    close(H)
    
    %%
    SingleIND = MAXMINZ(lists,1):MAXMINZ(lists,2);
    for i = SingleIND
        asize = size(DATSHOWDEFAULT,2);
        bsize = size(DATSHOWDEFAULT,1);        
        TempOut = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,i));
        Dsize = [size(TempOut,2),size(TempOut,1)];
        factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
        H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
        imagesc(TempOut,[0,dshownum*2]);colormap(colmaptemp);
        axis off;
        saveas(H,[OUTDIRSDIR1,filesep,'SepRegionSlice',num2str(i),'.fig'])
        set(H,'PaperPositionMode','manual');
        set(H,'PaperUnits','inch')
        XSIZE = size(TempOut,2);
        YSIZE = size(TempOut,1);
        factor = 1:100;
        XSIZEnew = XSIZE*factor;
        YSIZEnew = YSIZE*factor;
        FACTORS1 = find(XSIZEnew>Hsize(3));
        FACTORS2 = find(YSIZEnew>Hsize(4));
        FACTORS = max(FACTORS1(1),FACTORS2(1));
        XSIZEU = XSIZEnew(FACTORS);
        YSIZEU = YSIZEnew(FACTORS);
        set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
        print(H,[OUTDIRSDIR1,filesep,'SepRegionSlice',num2str(i),'.tif'],'-dtiff','-r300')
        close(H)
        for j = 1:dshownum
            Showtempusedtemp = rot90(DATSHOWDEFAULT(startpointY:endpointY,startpointX:endpointX,i));
            Showtempusedtemp2 = rot90(DATBG(startpointY:endpointY,startpointX:endpointX,i));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            Dsize = [size(Showtempusedtemp2,2),size(Showtempusedtemp2,1)];
            factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
            H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
            imagesc(Showtempusedtemp2,[0,dshownum*2]);colormap(colmaptemp);
            axis off;
            saveas(H,[OUTDIRSDIR1,filesep,'SepRegionSlice',num2str(i),'ROI',num2str(j),'.fig'])
            set(H,'PaperPositionMode','manual');
            set(H,'PaperUnits','inch')
            XSIZE = size(Showtempusedtemp2,2);
            YSIZE = size(Showtempusedtemp2,1);
            factor = 1:100;
            XSIZEnew = XSIZE*factor;
            YSIZEnew = YSIZE*factor;
            FACTORS1 = find(XSIZEnew>Hsize(3));
            FACTORS2 = find(YSIZEnew>Hsize(4));
            FACTORS = max(FACTORS1(1),FACTORS2(1));
            XSIZEU = XSIZEnew(FACTORS);
            YSIZEU = YSIZEnew(FACTORS);
            set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
            print(H,[OUTDIRSDIR1,filesep,'SepRegionSlice',num2str(i),'ROI',num2str(j),'.tif'],'-dtiff','-r300')
            close(H)
        end
    end
elseif corv
    asize = size(DATSHOWDEFAULT,3);
    bsize = size(DATSHOWDEFAULT,1);
    MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(squeeze(DATSHOWDEFAULT(:,sliceorderdefault(i),:)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(:,sliceorderdefault(i),:)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(:,sliceorderdefault(i),:)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
        end
    end
    Dsize = [size(MATSHOWMATRIX,2),size(MATSHOWMATRIX,1)];
    factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
    H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
    imagesc(MATSHOWMATRIX,[0,dshownum*2]);colormap(colmaptemp);
    axis off;
    saveas(H,[OUTDIRSDIR2,filesep,'CombinedWhole.fig'])
    set(H,'PaperPositionMode','manual');
    set(H,'PaperUnits','inch')
    XSIZE = size(MATSHOWMATRIX,2);
    YSIZE = size(MATSHOWMATRIX,1);
    factor = 1:100;
    XSIZEnew = XSIZE*factor;
    YSIZEnew = YSIZE*factor;
    FACTORS1 = find(XSIZEnew>Hsize(3));
    FACTORS2 = find(YSIZEnew>Hsize(4));
    FACTORS = max(FACTORS1(1),FACTORS2(1));
    XSIZEU = XSIZEnew(FACTORS);
    YSIZEU = YSIZEnew(FACTORS);
    set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
    print(H,[OUTDIRSDIR2,filesep,'CombinedWhole.tif'],'-dtiff','-r300')
    close(H)
    %%
    SingleIND = MAXMINY(lists,1):MAXMINY(lists,2);
    for i = SingleIND
        asize = size(DATSHOWDEFAULT,3);
        bsize = size(DATSHOWDEFAULT,1);        
        TempOut = rot90(squeeze(DATSHOWDEFAULT(:,i,:)));
        Dsize = [size(TempOut,2),size(TempOut,1)];
        factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
        H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
        imagesc(TempOut,[0,dshownum*2]);colormap(colmaptemp);
        axis off;
        saveas(H,[OUTDIRSDIR2,filesep,'SepWholeSlice',num2str(i),'.fig'])
        set(H,'PaperPositionMode','manual');
        set(H,'PaperUnits','inch')
        XSIZE = size(TempOut,2);
        YSIZE = size(TempOut,1);
        factor = 1:100;
        XSIZEnew = XSIZE*factor;
        YSIZEnew = YSIZE*factor;
        FACTORS1 = find(XSIZEnew>Hsize(3));
        FACTORS2 = find(YSIZEnew>Hsize(4));
        FACTORS = max(FACTORS1(1),FACTORS2(1));
        XSIZEU = XSIZEnew(FACTORS);
        YSIZEU = YSIZEnew(FACTORS);
        set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
        print(H,[OUTDIRSDIR2,filesep,'SepWholeSlice',num2str(i),'.tif'],'-dtiff','-r300')
        close(H)
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(:,i,:)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(:,i,:)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            Dsize = [size(Showtempusedtemp2,2),size(Showtempusedtemp2,1)];
            factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
            H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
            imagesc(Showtempusedtemp2,[0,dshownum*2]);colormap(colmaptemp);
            axis off;
            saveas(H,[OUTDIRSDIR2,filesep,'SepWholeSlice',num2str(i),'ROI',num2str(j),'.fig'])
            set(H,'PaperPositionMode','manual');
            set(H,'PaperUnits','inch')
            XSIZE = size(Showtempusedtemp2,2);
            YSIZE = size(Showtempusedtemp2,1);
            factor = 1:100;
            XSIZEnew = XSIZE*factor;
            YSIZEnew = YSIZE*factor;
            FACTORS1 = find(XSIZEnew>Hsize(3));
            FACTORS2 = find(YSIZEnew>Hsize(4));
            FACTORS = max(FACTORS1(1),FACTORS2(1));
            XSIZEU = XSIZEnew(FACTORS);
            YSIZEU = YSIZEnew(FACTORS);
            set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
            print(H,[OUTDIRSDIR2,filesep,'SepWholeSlice',num2str(i),'ROI',num2str(j),'.tif'],'-dtiff','-r300')
            close(H)
        end
    end
    %%
    POSX = MAXMINZ(lists,:);
    POSY = MAXMINX(lists,:);
    LENX = abs(POSX(1)-POSX(2));
    LENY = abs(POSY(1)-POSY(2));
    FACTTEMP = get(ASSHOWMAP.factoredit,'string');
    factor = str2num(FACTTEMP); % This is the change for the cut show.
    FACTOR = 0.1;
%     ASSHOWMAP.factor = factor;
%     ASSHOWMAP.FACTOR = FACTOR;
    sepfactorX = ceil(LENX*FACTOR);
    sepfactorY = ceil(LENY*FACTOR);
    halffactor = factor/2;
    LENXh = ceil(LENX*halffactor);
    LENYh = ceil(LENY*halffactor);
    startpointX = POSX(1)-LENXh;
    startpointY = POSY(1)-LENYh;
    endpointX = POSX(2)+LENXh;
    endpointY = POSY(2)+LENYh;
    if startpointY<1
        startpointY = 1;
    end
    if startpointX<1
        startpointX = 1;
    end
    if endpointX>asize
        endpointX = asize;
    end
    if endpointY>bsize
        endpointY = bsize;
    end
    asizenew = length(startpointX:endpointX);
    bsizenew = length(startpointY:endpointY);
    MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(squeeze(DATSHOWDEFAULT(startpointY:endpointY,sliceorderdefault(i),startpointX:endpointX)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(startpointY:endpointY,sliceorderdefault(i),startpointX:endpointX)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(startpointY:endpointY,sliceorderdefault(i),startpointX:endpointX)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX2(1+asizenew*j+sepfactorX*(j+1):asizenew*(j+1)+sepfactorX*(j+1),1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
        end
    end
    Dsize = [size(MATSHOWMATRIX2,2),size(MATSHOWMATRIX2,1)];
    factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
    H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
    imagesc(MATSHOWMATRIX2,[0,dshownum*2]);colormap(colmaptemp);
    axis off;
    saveas(H,[OUTDIRSDIR2,filesep,'CombinedRegion.fig'])
    set(H,'PaperPositionMode','manual');
    set(H,'PaperUnits','inch')
    XSIZE = size(MATSHOWMATRIX2,2);
    YSIZE = size(MATSHOWMATRIX2,1);
    factor = 1:100;
    XSIZEnew = XSIZE*factor;
    YSIZEnew = YSIZE*factor;
    FACTORS1 = find(XSIZEnew>Hsize(3));
    FACTORS2 = find(YSIZEnew>Hsize(4));
    FACTORS = max(FACTORS1(1),FACTORS2(1));
    XSIZEU = XSIZEnew(FACTORS);
    YSIZEU = YSIZEnew(FACTORS);
    set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
    print(H,[OUTDIRSDIR2,filesep,'CombinedRegion.tif'],'-dtiff','-r300')
    close(H)    
    %%
    SingleIND = MAXMINY(lists,1):MAXMINY(lists,2);
    for i = SingleIND
        asize = size(DATSHOWDEFAULT,2);
        bsize = size(DATSHOWDEFAULT,1);        
        TempOut = rot90(squeeze(DATSHOWDEFAULT(startpointY:endpointY,i,startpointX:endpointX)));
        Dsize = [size(TempOut,2),size(TempOut,1)];
        factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
        H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
        imagesc(TempOut,[0,dshownum*2]);colormap(colmaptemp);
        axis off;
        saveas(H,[OUTDIRSDIR2,filesep,'SepRegionSlice',num2str(i),'.fig'])
        set(H,'PaperPositionMode','manual');
        set(H,'PaperUnits','inch')
        XSIZE = size(TempOut,2);
        YSIZE = size(TempOut,1);
        factor = 1:100;
        XSIZEnew = XSIZE*factor;
        YSIZEnew = YSIZE*factor;
        FACTORS1 = find(XSIZEnew>Hsize(3));
        FACTORS2 = find(YSIZEnew>Hsize(4));
        FACTORS = max(FACTORS1(1),FACTORS2(1));
        XSIZEU = XSIZEnew(FACTORS);
        YSIZEU = YSIZEnew(FACTORS);
        set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
        print(H,[OUTDIRSDIR2,filesep,'SepRegionSlice',num2str(i),'.tif'],'-dtiff','-r300')
        close(H)
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(startpointY:endpointY,i,startpointX:endpointX)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(startpointY:endpointY,i,startpointX:endpointX)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            Dsize = [size(Showtempusedtemp2,2),size(Showtempusedtemp2,1)];
            factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
            H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
            imagesc(Showtempusedtemp2,[0,dshownum*2]);colormap(colmaptemp);
            axis off;
            saveas(H,[OUTDIRSDIR2,filesep,'SepRegionSlice',num2str(i),'ROI',num2str(j),'.fig'])
            set(H,'PaperPositionMode','manual');
            set(H,'PaperUnits','inch')
            XSIZE = size(Showtempusedtemp2,2);
            YSIZE = size(Showtempusedtemp2,1);
            factor = 1:100;
            XSIZEnew = XSIZE*factor;
            YSIZEnew = YSIZE*factor;
            FACTORS1 = find(XSIZEnew>Hsize(3));
            FACTORS2 = find(YSIZEnew>Hsize(4));
            FACTORS = max(FACTORS1(1),FACTORS2(1));
            XSIZEU = XSIZEnew(FACTORS);
            YSIZEU = YSIZEnew(FACTORS);
            set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
            print(H,[OUTDIRSDIR2,filesep,'SepRegionSlice',num2str(i),'ROI',num2str(j),'.tif'],'-dtiff','-r300')
            close(H)
        end
    end
elseif sagv
    asize = size(DATSHOWDEFAULT,3);
    bsize = size(DATSHOWDEFAULT,2);
    MATSHOWMATRIX = zeros(asize*(dshownum+1),bsize*length(sliceorderdefault));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX(1:asize,1+bsize*(i-1):bsize*i) = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),:,:)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),:,:)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(sliceorderdefault(i),:,:)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX(1+asize*j:asize*(j+1),1+bsize*(i-1):bsize*i) = Showtempusedtemp2;
        end
    end
    Dsize = [size(MATSHOWMATRIX,2),size(MATSHOWMATRIX,1)];
    factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
    H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
    imagesc(MATSHOWMATRIX,[0,dshownum*2]);colormap(colmaptemp);
    axis off;
    saveas(H,[OUTDIRSDIR3,filesep,'CombinedWhole.fig'])
    set(H,'PaperPositionMode','manual');
    set(H,'PaperUnits','inch')
    XSIZE = size(MATSHOWMATRIX,2);
    YSIZE = size(MATSHOWMATRIX,1);
    factor = 1:100;
    XSIZEnew = XSIZE*factor;
    YSIZEnew = YSIZE*factor;
    FACTORS1 = find(XSIZEnew>Hsize(3));
    FACTORS2 = find(YSIZEnew>Hsize(4));
    FACTORS = max(FACTORS1(1),FACTORS2(1));
    XSIZEU = XSIZEnew(FACTORS);
    YSIZEU = YSIZEnew(FACTORS);
    set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
    print(H,[OUTDIRSDIR3,filesep,'CombinedWhole.tif'],'-dtiff','-r300')
    close(H)
    %%
    SingleIND = MAXMINX(lists,1):MAXMINX(lists,2);
    for i = SingleIND
        asize = size(DATSHOWDEFAULT,3);
        bsize = size(DATSHOWDEFAULT,2);        
        TempOut = rot90(squeeze(DATSHOWDEFAULT(i,:,:)));
        Dsize = [size(TempOut,2),size(TempOut,1)];
        factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
        H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
        imagesc(TempOut,[0,dshownum*2]);colormap(colmaptemp);
        axis off;
        saveas(H,[OUTDIRSDIR3,filesep,'SepWholeSlice',num2str(i),'.fig'])
        set(H,'PaperPositionMode','manual');
        set(H,'PaperUnits','inch')
        XSIZE = size(TempOut,2);
        YSIZE = size(TempOut,1);
        factor = 1:100;
        XSIZEnew = XSIZE*factor;
        YSIZEnew = YSIZE*factor;
        FACTORS1 = find(XSIZEnew>Hsize(3));
        FACTORS2 = find(YSIZEnew>Hsize(4));
        FACTORS = max(FACTORS1(1),FACTORS2(1));
        XSIZEU = XSIZEnew(FACTORS);
        YSIZEU = YSIZEnew(FACTORS);
        set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
        print(H,[OUTDIRSDIR3,filesep,'SepWholeSlice',num2str(i),'.tif'],'-dtiff','-r300')
        close(H)
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(i,:,:)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(i,:,:)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            Dsize = [size(Showtempusedtemp2,2),size(Showtempusedtemp2,1)];
            factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
            H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
            imagesc(Showtempusedtemp2,[0,dshownum*2]);colormap(colmaptemp);
            axis off;
            saveas(H,[OUTDIRSDIR3,filesep,'SepWholeSlice',num2str(i),'ROI',num2str(j),'.fig'])
            set(H,'PaperPositionMode','manual');
            set(H,'PaperUnits','inch')
            XSIZE = size(Showtempusedtemp2,2);
            YSIZE = size(Showtempusedtemp2,1);
            factor = 1:100;
            XSIZEnew = XSIZE*factor;
            YSIZEnew = YSIZE*factor;
            FACTORS1 = find(XSIZEnew>Hsize(3));
            FACTORS2 = find(YSIZEnew>Hsize(4));
            FACTORS = max(FACTORS1(1),FACTORS2(1));
            XSIZEU = XSIZEnew(FACTORS);
            YSIZEU = YSIZEnew(FACTORS);
            set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
            print(H,[OUTDIRSDIR3,filesep,'SepWholeSlice',num2str(i),'ROI',num2str(j),'.tif'],'-dtiff','-r300')
            close(H)
        end
    end
    %%
    POSX = MAXMINZ(lists,:);
    POSY = MAXMINY(lists,:);
    LENX = abs(POSX(1)-POSX(2));
    LENY = abs(POSY(1)-POSY(2));
    FACTTEMP = get(ASSHOWMAP.factoredit,'string');
    factor = str2num(FACTTEMP); % This is the change for the cut show.
    FACTOR = 0.1;
    ASSHOWMAP.factor = factor;
    ASSHOWMAP.FACTOR = FACTOR;
    sepfactorX = ceil(LENX*FACTOR);
    sepfactorY = ceil(LENY*FACTOR);
    halffactor = factor/2;
    LENXh = ceil(LENX*halffactor);
    LENYh = ceil(LENY*halffactor);
    startpointX = POSX(1)-LENXh;
    startpointY = POSY(1)-LENYh;
    endpointX = POSX(2)+LENXh;
    endpointY = POSY(2)+LENYh;
    if startpointY<1
        startpointY = 1;
    end
    if startpointX<1
        startpointX = 1;
    end
    if endpointX>asize
        endpointX = asize;
    end
    if endpointY>bsize
        endpointY = bsize;
    end
    asizenew = length(startpointX:endpointX);
    bsizenew = length(startpointY:endpointY);
    MATSHOWMATRIX2 = zeros(asizenew*(dshownum+1)+sepfactorX*(dshownum+2),bsizenew*length(sliceorderdefault)+sepfactorY*(length(sliceorderdefault)+1));
    for i = 1:length(sliceorderdefault)
        MATSHOWMATRIX2(1+sepfactorX:asizenew+sepfactorX,1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),startpointY:endpointY,startpointX:endpointX)));
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(sliceorderdefault(i),startpointY:endpointY,startpointX:endpointX)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(sliceorderdefault(i),startpointY:endpointY,startpointX:endpointX)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            MATSHOWMATRIX2(1+asizenew*j+sepfactorX*(j+1):asizenew*(j+1)+sepfactorX*(j+1),1+sepfactorY*i+bsizenew*(i-1):bsizenew*i+sepfactorY*i) = Showtempusedtemp2;
        end
    end
    
    Dsize = [size(MATSHOWMATRIX2,2),size(MATSHOWMATRIX2,1)];
    factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
    H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
    imagesc(MATSHOWMATRIX2,[0,dshownum*2]);colormap(colmaptemp);
    axis off;
    saveas(H,[OUTDIRSDIR3,filesep,'CombinedRegion.fig'])
    set(H,'PaperPositionMode','manual');
    set(H,'PaperUnits','inch')
    XSIZE = size(MATSHOWMATRIX2,2);
    YSIZE = size(MATSHOWMATRIX2,1);
    factor = 1:100;
    XSIZEnew = XSIZE*factor;
    YSIZEnew = YSIZE*factor;
    FACTORS1 = find(XSIZEnew>Hsize(3));
    FACTORS2 = find(YSIZEnew>Hsize(4));
    FACTORS = max(FACTORS1(1),FACTORS2(1));
    XSIZEU = XSIZEnew(FACTORS);
    YSIZEU = YSIZEnew(FACTORS);
    set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
    print(H,[OUTDIRSDIR3,filesep,'CombinedRegion.tif'],'-dtiff','-r300')
    close(H)
    %%
    SingleIND = MAXMINX(lists,1):MAXMINX(lists,2);
    for i = SingleIND
        asize = size(DATSHOWDEFAULT,2);
        bsize = size(DATSHOWDEFAULT,1);        
        TempOut = rot90(squeeze(DATSHOWDEFAULT(i,startpointY:endpointY,startpointX:endpointX)));
        Dsize = [size(TempOut,2),size(TempOut,1)];
        factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
        H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
        imagesc(TempOut,[0,dshownum*2]);colormap(colmaptemp);
        axis off;
        saveas(H,[OUTDIRSDIR3,filesep,'SepRegionSlice',num2str(i),'.fig'])
        set(H,'PaperPositionMode','manual');
        set(H,'PaperUnits','inch')
        XSIZE = size(TempOut,2);
        YSIZE = size(TempOut,1);
        factor = 1:100;
        XSIZEnew = XSIZE*factor;
        YSIZEnew = YSIZE*factor;
        FACTORS1 = find(XSIZEnew>Hsize(3));
        FACTORS2 = find(YSIZEnew>Hsize(4));
        FACTORS = max(FACTORS1(1),FACTORS2(1));
        XSIZEU = XSIZEnew(FACTORS);
        YSIZEU = YSIZEnew(FACTORS);
        set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
        print(H,[OUTDIRSDIR3,filesep,'SepRegionSlice',num2str(i),'.tif'],'-dtiff','-r300')
        close(H)
        for j = 1:dshownum
            Showtempusedtemp = rot90(squeeze(DATSHOWDEFAULT(i,startpointY:endpointY,startpointX:endpointX)));
            Showtempusedtemp2 = rot90(squeeze(DATBG(i,startpointY:endpointY,startpointX:endpointX)));
            Showtempusedtemp2(Showtempusedtemp==dshownum-0.5+dmidindex(j+1)) = dshownum-0.5+dmidindex(j+1);
            Dsize = [size(Showtempusedtemp2,2),size(Showtempusedtemp2,1)];
            factorD = min(Dsize(1)/Hexist2,Dsize(2)/Hexist1);
            H = figure('pos',[100,100,Dsize(1)*factorD,Dsize(2)*factorD]);
            imagesc(Showtempusedtemp2,[0,dshownum*2]);colormap(colmaptemp);
            axis off;
            saveas(H,[OUTDIRSDIR3,filesep,'SepRegionSlice',num2str(i),'ROI',num2str(j),'.fig'])
            set(H,'PaperPositionMode','manual');
            set(H,'PaperUnits','inch')
            XSIZE = size(Showtempusedtemp2,2);
            YSIZE = size(Showtempusedtemp2,1);
            factor = 1:100;
            XSIZEnew = XSIZE*factor;
            YSIZEnew = YSIZE*factor;
            FACTORS1 = find(XSIZEnew>Hsize(3));
            FACTORS2 = find(YSIZEnew>Hsize(4));
            FACTORS = max(FACTORS1(1),FACTORS2(1));
            XSIZEU = XSIZEnew(FACTORS);
            YSIZEU = YSIZEnew(FACTORS);
            set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
            print(H,[OUTDIRSDIR3,filesep,'SepRegionSlice',num2str(i),'ROI',num2str(j),'.tif'],'-dtiff','-r300')
            close(H)
        end
    end
end
end
function Matrixmap(varargin)
ASSHOWMAP = varargin{3};
pat = ASSHOWMAP.pat;
load(fullfile(pat,'SEEDCOLOR.mat'));
load(fullfile(pat,'graycolmap.mat'));
targetnum = ASSHOWMAP.Info.targetnum;
targetdir = ASSHOWMAP.Info.target;
dshownum = ASSHOWMAP.dshownum;
% dindir = fullfile(ASSHOWMAP.Info.input,'mixedMaxID.nii');
% [vmix dmix] = Dynamic_read_dir_NIFTI(dindir);
[vtar dtar] = Dynamic_read_dir_NIFTI(targetdir);
ind = find(dtar);
load(fullfile(ASSHOWMAP.Info.input,'computeval.mat'));
figure('units','norm','pos',[0.2 0.45 0.7 0.4]);
imagesc(r,[-0.5 0.5]);colormap('jet')
title('Original Rmap')
[maxv maxid] = max(r);
figure('units','norm','pos',[0.2 0.15 0.7 0.2]);
imagesc(maxid,[0 dshownum+0.5]);colormap(SEEDCOLORSHOW(1:dshownum,:))
title('Original WTA value')
[ix iy] = sort(maxid);
figure('units','norm','pos',[0.2 0.45 0.7 0.4]);
imagesc(r(:,iy),[-0.5 0.5]);colormap('jet');
title('Sorted Rmap')
figure('units','norm','pos',[0.2 0.15 0.7 0.2]);
imagesc(ix,[0 dshownum+0.5]);colormap(SEEDCOLORSHOW(1:dshownum,:))
title('Sorted WTA value')
DTARU = dtar(ind);
dtaruind = unique(DTARU);
for i = 1:targetnum
    INDTEMP = (DTARU==dtaruind(i));
    Rtemp = r(:,INDTEMP);
    figure('units','norm','pos',[0.2 0.45 0.7 0.4]);
    imagesc(Rtemp,[-0.5 0.5]);colormap('jet')
    title(['ROI',num2str(i),' Rmap'])
    [maxv maxid] = max(Rtemp);
    figure('units','norm','pos',[0.2 0.15 0.7 0.2]);
    imagesc(maxid,[0 dshownum+0.5]);colormap(SEEDCOLORSHOW(1:dshownum,:))
    title(['ROI',num2str(i),' WTA value'])
    [ix iy] = sort(maxid);
    figure('units','norm','pos',[0.2 0.45 0.7 0.4]);
    imagesc(Rtemp(:,iy),[-0.5 0.5]);colormap('jet');
    title(['ROI',num2str(i),' Sorted Rmap'])
    figure('units','norm','pos',[0.2 0.15 0.7 0.2]);
    imagesc(ix,[0 dshownum+0.5]);colormap(SEEDCOLORSHOW(1:dshownum,:))
    title(['ROI',num2str(i),' Sorted WTA value'])
end

end

%%
function CirclePlot(R,ax)
alpha=0:pi/20:2*pi;    %�Ƕ�[0,2*pi]

% R=2;                   %�뾶
x=R*cos(alpha);
y=R*sin(alpha);
plot(x,y,'--','parent',ax,'linewidth',0.5,'color',[0.5,0.5,0.5]);
end
function Outv = FindMaxNearVal(Inv)
valtry = [0.000001 0.00001 0.0001 0.001 0.01 0.1 1 10 100 1000 10000 100000 1000000];
ind = ones(7,1);
for i = 1:length(valtry)
    valv = Inv/valtry(i);
    if valv>1
        ind(i) = 0;
    end
end
indmin = find(ind==0);
indsel = indmin(end);
VALTRY = valtry(indsel);
Outv = ceil((Inv/VALTRY)*10)/10*VALTRY;
end
function Outv = FindMaxNearValper(Inv)
valtry = [0.000001 0.00001 0.0001 0.001 0.01 0.1 1 10 100 1000 10000 100000 1000000];
ind = ones(7,1);
for i = 1:length(valtry)
    valv = Inv/valtry(i);
    if valv>1
        ind(i) = 0;
    end
end
indmin = find(ind==0);
indsel = indmin(end);
VALTRY = valtry(indsel+1);
Outv = ceil((Inv/VALTRY)*10)/10*VALTRY;
end
function Radbg(R,ax,dshownum)
for i = 1:dshownum
    AlphaVal(i,1) = pi/2-2*pi/dshownum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)]*R,[0,YShow(i,1)]*R,'k--','parent',ax);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
for i = 1:dshownum
    for j = 1:5
        plot([XShow2(i,1)*j/5,XShow2(i+1,1)*j/5]*R,...
            [YShow2(i,1)*j/5,YShow2(i+1,1)*j/5]*R,...
            '--','parent',ax,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
end
if R>=5
    for j = 1:5
        text(XShow2(1,1)*j/5*R,YShow2(1,1)*j/5*R,num2str(round(R/5*j)),'parent',ax)
    end
else    
    for j = 1:5
        text(XShow2(1,1)*j/5*R,YShow2(1,1)*j/5*R,num2str(R/5*j),'parent',ax)
    end
end
    
end
function Radbg2(R,ax,dshownum)
for i = 1:dshownum
    AlphaVal(i,1) = pi/2-2*pi/dshownum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
    plot([0,XShow(i,1)]*R,[0,YShow(i,1)]*R,'k--','parent',ax);
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
for i = 1:dshownum
    for j = 1:5
        plot([XShow2(i,1)*j/5,XShow2(i+1,1)*j/5]*R,...
            [YShow2(i,1)*j/5,YShow2(i+1,1)*j/5]*R,...
            '--','parent',ax,'linewidth',0.5,'color',[0.5,0.5,0.5]);
    end
end
% for j = 1:5
%     text(XShow2(1,1)*j/5*R,YShow2(1,1)*j/5*R,num2str(round(R/5*j)),'parent',ax)
% end    
end
function PlotRadar(Patternnum,Hax,SEEDCOLORSHOW,dshownum)
if size(Patternnum,1)>size(Patternnum,2)
    Patternnum = Patternnum';
end
Patternnum2 = [Patternnum,Patternnum(1)];
for i = 1:dshownum
    AlphaVal(i,1) = pi/2-2*pi/dshownum*(i-1);
    XShow(i,1) = cos(AlphaVal(i,1));
    YShow(i,1) = sin(AlphaVal(i,1));
end
XShow2 = [XShow;XShow(1)];
YShow2 = [YShow;YShow(1)];
for i = 1:dshownum
    plot([XShow2(i,1)*Patternnum2(i),XShow2(i+1,1)*Patternnum2(i+1)],...
        [YShow2(i,1)*Patternnum2(i),YShow2(i+1,1)*Patternnum2(i+1)],...
        'parent',Hax,...
        'linewidth',2,...
        'color','k');
end
for i = 1:dshownum   
    plot(XShow2(i,1)*Patternnum2(i),YShow2(i,1)*Patternnum2(i),'parent',Hax,...
        'marker','o',...
        'markersize',8,'MarkerFaceColor',SEEDCOLORSHOW(i,:),...
        'MarkerEdgeColor',SEEDCOLORSHOW(i,:)) 
end
end