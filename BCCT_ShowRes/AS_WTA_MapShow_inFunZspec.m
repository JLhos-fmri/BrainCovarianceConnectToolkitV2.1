function AS_WTA_MapShow_inFunZspec(inmixdir,outdir,width,lenth)
[pat nam ext] = fileparts(which('AS_WTA_MapShow.m'));
backgroundmap = fullfile(pat,'mni_icbm152_t1_tal_nlin_asym_09a.nii');
% indir = uigetdir(pwd,'Result of WTA');
% outdir = uigetdir(pwd,'Pic Out Dir');
[vbg dbg] = Dynamic_read_dir_NIFTI(backgroundmap);
dbgre = reshape(dbg,vbg.dim(1),vbg.dim(2),vbg.dim(3));
% inmixdir = fullfile(indir,'mixedMaxID.nii');
[vmid dmid] = Dynamic_read_dir_NIFTI(inmixdir);
load(fullfile(pat,'SEEDCOLOR.mat'));
load(fullfile(pat,'graycolmap.mat'));
dmidindex = unique(dmid);
dshownum = length(dmidindex)-1;
dmidre = reshape(dmid,vmid.dim(1),vmid.dim(2),vmid.dim(3));
indnew = find(dmidre>0);
[indnewx, indnewy, indnewz] = ind2sub(vmid.dim,indnew);

indnewsub = [indnewx, indnewy, indnewz];
indnewsubT = trans2solution(indnewsub,vmid.mat,vbg.mat);
diffx = abs(vmid.mat(1,1)/vbg.mat(1,1))/2*abs(vbg.mat(1,1));
diffy = abs(vmid.mat(2,2)/vbg.mat(2,2))/2*abs(vbg.mat(2,2));
diffz = abs(vmid.mat(3,3)/vbg.mat(3,3))/2*abs(vbg.mat(3,3));
% INDNEWPLOT = [];
MATOUTT = zeros(vbg.dim);
for i = 1:size(indnewsub,1)
    indu = indnewsubT(i,:);
    induextendx = round(indu(1)-diffx):round(indu(1)+diffx);
    induextendy = round(indu(2)-diffy):round(indu(2)+diffy);
    induextendz = round(indu(3)-diffz):round(indu(3)+diffz);
    MATOUTT(induextendx,induextendy,induextendz) = dmidre(indnew(i));
end

INDSHOW = find(MATOUTT>0);
[ix iy iz] = ind2sub(vbg.dim,INDSHOW);

rangx = unique(ix);
rangy = unique(iy);
rangz = unique(iz);
% mnixyz = cor2mni([ix,iy iz],vbg.mat);
xranglen = abs(ix(end)-ix(1));
yranglen = abs(iy(end)-iy(1));
zranglen = abs(iz(end)-iz(1));
indsep = round(1:63/dshownum:64);
colmaptemp = [];
for i = 1:dshownum
    colmaptemp(indsep(i):indsep(i+1),:) = ones(length(indsep(i):indsep(i+1)),1)*[SEEDCOLORSHOW(i,1),SEEDCOLORSHOW(i,2),SEEDCOLORSHOW(i,3)];
end
colmaptemp = [colgray;colmaptemp];
% for i = 1:10
for i = 1:length(rangz)
    Doutshowtemp = squeeze(dbgre(:,:,rangz(i)));
    Doutshowtemp(Doutshowtemp<30*0.9) = 0;
    Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
    Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
    Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum);
%     figure;imagesc(Doutshowtemp)
    Doutmatout = squeeze(MATOUTT(:,:,rangz(i)));
    for j = 1:dshownum
        ind = find(Doutmatout==dmidindex(j+1));
        Doutshowtemp(ind) = dshownum+0.5+(j-1);
    end
    Doutshowtemp1 = rot90(Doutshowtemp);
    H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
    imagesc(Doutshowtemp1,[0,dshownum*2]);colormap(colmaptemp);
    axis off;
    saveas(H,[outdir,filesep,['z',num2str(rangz(i)),'.fig']])
    set(H,'PaperPositionMode','manual');
    set(H,'PaperUnits','inch')
    XSIZE = size(Doutshowtemp1,2);
    YSIZE = size(Doutshowtemp1,1);
    factor = 1:100;
    XSIZEnew = XSIZE*factor;
    YSIZEnew = YSIZE*factor;
    FACTORS1 = find(XSIZEnew>400);
    FACTORS2 = find(YSIZEnew>400);
    FACTORS = max(FACTORS1(1),FACTORS2(1));
    XSIZEU = XSIZEnew(FACTORS);
    YSIZEU = YSIZEnew(FACTORS);    
    set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
    print(H,[outdir,filesep,['z',num2str(rangz(i)),'.tif']],'-dtiff','-r300')
    close(H)
end

% for i = 1:length(rangy)
%     Doutshowtemp = squeeze(dbgre(:,rangy(i),:));
%     Doutshowtemp(Doutshowtemp<30*0.9) = 0;
%     Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
%     Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
%     Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum);
% %     figure;imagesc(Doutshowtemp)
%     Doutmatout = squeeze(MATOUTT(:,rangy(i),:));
%     for j = 1:dshownum
%         ind = find(Doutmatout==dmidindex(j+1));
%         Doutshowtemp(ind) = dshownum+0.5+(j-1);
%     end
%     Doutshowtemp1 = rot90(Doutshowtemp);
%     H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
%     imagesc(Doutshowtemp1,[0,dshownum*2]);colormap(colmaptemp);
%     axis off;
%     saveas(H,[outdir,filesep,['y',num2str(rangy(i)),'.fig']])
%     set(H,'PaperPositionMode','manual');
%     set(H,'PaperUnits','inch')
%     XSIZE = size(Doutshowtemp1,2);
%     YSIZE = size(Doutshowtemp1,1);
%     factor = 1:100;
%     XSIZEnew = XSIZE*factor;
%     YSIZEnew = YSIZE*factor;
%     FACTORS1 = find(XSIZEnew>400);
%     FACTORS2 = find(YSIZEnew>400);
%     FACTORS = max(FACTORS1(1),FACTORS2(1));
%     XSIZEU = XSIZEnew(FACTORS);
%     YSIZEU = YSIZEnew(FACTORS);    
%     set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
%     print(H,[outdir,filesep,['y',num2str(rangy(i)),'.tif']],'-dtiff','-r300')
%     close(H)
% end

% for i = 1:length(rangx)
%     Doutshowtemp = squeeze(dbgre(rangx(i),:,:));
%     Doutshowtemp(Doutshowtemp<30*0.9) = 0;
%     Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
%     Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
%     Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum);
% %     figure;imagesc(Doutshowtemp)
%     Doutmatout = squeeze(MATOUTT(rangx(i),:,:));
%     for j = 1:dshownum
%         ind = find(Doutmatout==dmidindex(j+1));
%         Doutshowtemp(ind) = dshownum+0.5+(j-1);
%     end
%     Doutshowtemp1 = rot90(Doutshowtemp);
%     H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
%     imagesc(Doutshowtemp1,[0,dshownum*2]);colormap(colmaptemp);
%     axis off;
%     saveas(H,[outdir,filesep,['x',num2str(rangx(i)),'.fig']])
%     set(H,'PaperPositionMode','manual');
%     set(H,'PaperUnits','inch')
%     XSIZE = size(Doutshowtemp1,2);
%     YSIZE = size(Doutshowtemp1,1);
%     factor = 1:100;
%     XSIZEnew = XSIZE*factor;
%     YSIZEnew = YSIZE*factor;
%     FACTORS1 = find(XSIZEnew>400);
%     FACTORS2 = find(YSIZEnew>400);
%     FACTORS = max(FACTORS1(1),FACTORS2(1));
%     XSIZEU = XSIZEnew(FACTORS);
%     YSIZEU = YSIZEnew(FACTORS);    
%     set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
%     print(H,[outdir,filesep,['x',num2str(rangx(i)),'.tif']],'-dtiff','-r300')
%     close(H)
% end


%%
minx = floor((max(rangx)+min(rangx))/2-width/2);
maxx = ceil((max(rangx)+min(rangx))/2+width/2);
miny = floor((max(rangy)+min(rangy))/2-lenth/2);
maxy = ceil((max(rangy)+min(rangy))/2+lenth/2);
% minz = min(rangz)-(max(rangz)-min(rangz))*extendfactor3;
% maxz = max(rangz)+(max(rangz)-min(rangz))*extendfactor3;
for i = 1:length(rangz)
    Doutshowtemp = squeeze(dbgre(minx:maxx,miny:maxy,rangz(i)));
    Doutshowtemp(Doutshowtemp<30*0.9) = 0;
    Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
    Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
    Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum);
%     figure;imagesc(Doutshowtemp)
    Doutmatout = squeeze(MATOUTT(minx:maxx,miny:maxy,rangz(i)));
    for j = 1:dshownum
        ind = find(Doutmatout==dmidindex(j+1));
        Doutshowtemp(ind) = dshownum+0.5+(j-1);
    end
    Doutshowtemp1 = rot90(Doutshowtemp);
    H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
    imagesc(Doutshowtemp1,[0,dshownum*2]);colormap(colmaptemp);
    axis off;
    saveas(H,[outdir,filesep,['Sz',num2str(rangz(i)),'.fig']])
    set(H,'PaperPositionMode','manual');
    set(H,'PaperUnits','inch')
    XSIZE = size(Doutshowtemp1,2);
    YSIZE = size(Doutshowtemp1,1);
    factor = 1:100;
    XSIZEnew = XSIZE*factor;
    YSIZEnew = YSIZE*factor;
    FACTORS1 = find(XSIZEnew>400);
    FACTORS2 = find(YSIZEnew>400);
    FACTORS = max(FACTORS1(1),FACTORS2(1));
    XSIZEU = XSIZEnew(FACTORS);
    YSIZEU = YSIZEnew(FACTORS);    
    set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
    print(H,[outdir,filesep,['Sz',num2str(rangz(i)),'.tif']],'-dtiff','-r300')
    close(H)
end
% minx:maxx,miny:maxy,minz:maxz;
% for i = 1:length(rangy)
%     Doutshowtemp = squeeze(dbgre(minx:maxx,rangy(i),minz:maxz));
%     Doutshowtemp(Doutshowtemp<30*0.9) = 0;
%     Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
%     Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
%     Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum);
% %     figure;imagesc(Doutshowtemp)
%     Doutmatout = squeeze(MATOUTT(minx:maxx,rangy(i),minz:maxz));
%     for j = 1:dshownum
%         ind = find(Doutmatout==dmidindex(j+1));
%         Doutshowtemp(ind) = dshownum+0.5+(j-1);
%     end
%     Doutshowtemp1 = rot90(Doutshowtemp);
%     H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
%     imagesc(Doutshowtemp1,[0,dshownum*2]);colormap(colmaptemp);
%     axis off;
%     saveas(H,[outdir,filesep,['Sy',num2str(rangy(i)),'.fig']])
%     set(H,'PaperPositionMode','manual');
%     set(H,'PaperUnits','inch')
%     XSIZE = size(Doutshowtemp1,2);
%     YSIZE = size(Doutshowtemp1,1);
%     factor = 1:100;
%     XSIZEnew = XSIZE*factor;
%     YSIZEnew = YSIZE*factor;
%     FACTORS1 = find(XSIZEnew>400);
%     FACTORS2 = find(YSIZEnew>400);
%     FACTORS = max(FACTORS1(1),FACTORS2(1));
%     XSIZEU = XSIZEnew(FACTORS);
%     YSIZEU = YSIZEnew(FACTORS);    
%     set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
%     print(H,[outdir,filesep,['Sy',num2str(rangy(i)),'.tif']],'-dtiff','-r300')
%     close(H)
% end

% minx:maxx,miny:maxy,minz:maxz;
% for i = 1:length(rangx)
%     Doutshowtemp = squeeze(dbgre(rangx(i),miny:maxy,minz:maxz));
%     Doutshowtemp(Doutshowtemp<30*0.9) = 0;
%     Doutshowtemp(Doutshowtemp>0&Doutshowtemp<=30*0.9) = 30+(dshownum)/(85-30);
%     Doutshowtemp(Doutshowtemp>=85) = 85-dshownum/55;
%     Doutshowtemp(Doutshowtemp>0) = (Doutshowtemp(Doutshowtemp>0)-30)/(85-30)*(dshownum);
% %     figure;imagesc(Doutshowtemp)
%     Doutmatout = squeeze(MATOUTT(rangx(i),miny:maxy,minz:maxz));
%     for j = 1:dshownum
%         ind = find(Doutmatout==dmidindex(j+1));
%         Doutshowtemp(ind) = dshownum+0.5+(j-1);
%     end
%     Doutshowtemp1 = rot90(Doutshowtemp);
%     H = figure('pos',[100,100,size(Doutshowtemp1,2)*2,size(Doutshowtemp1,1)*2]);
%     imagesc(Doutshowtemp1,[0,dshownum*2]);colormap(colmaptemp);
%     axis off;
%     saveas(H,[outdir,filesep,['Sx',num2str(rangx(i)),'.fig']])
%     set(H,'PaperPositionMode','manual');
%     set(H,'PaperUnits','inch')
%     XSIZE = size(Doutshowtemp1,2);
%     YSIZE = size(Doutshowtemp1,1);
%     factor = 1:100;
%     XSIZEnew = XSIZE*factor;
%     YSIZEnew = YSIZE*factor;
%     FACTORS1 = find(XSIZEnew>400);
%     FACTORS2 = find(YSIZEnew>400);
%     FACTORS = max(FACTORS1(1),FACTORS2(1));
%     XSIZEU = XSIZEnew(FACTORS);
%     YSIZEU = YSIZEnew(FACTORS);    
%     set(H,'Paperposition',[1,1,XSIZEU/300,YSIZEU/300]);
%     print(H,[outdir,filesep,['Sx',num2str(rangx(i)),'.tif']],'-dtiff','-r300')
%     close(H)
% end

% RCP = load(fullfile(indir,'RealCompute.mat'));
% TargetDir = RCP.RealCompPara.Targetdir;
% [vtar dtar] = Dynamic_read_dir_NIFTI(TargetDir);
% 
% tarnum = unique(dtar);
% for i = 
