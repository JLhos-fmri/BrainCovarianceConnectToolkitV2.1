function Show_GC_res_subfun(R,P3,NROI,colormapshow,enhanind,outdir,outname,width1,width2)
Hsize = get(0,'ScreenSize');
Bsize = min(Hsize(3),Hsize(4))*0.9;
Bsize2 = Bsize;
POSm1 = [20,20,Bsize2,Bsize2];
Hm1 = figure('pos',POSm1);
maxv = max(R(:));
axHm1_1 = axes('parent',Hm1,'unit','norm','pos',[0.1 0.1 0.8 0.8],'CLim',[-maxv maxv]);
image(R.*P3,'parent',axHm1_1,'CDataMapping','scaled');
axis(axHm1_1,[-0.5,NROI+1.5,-0.5,NROI+1.5])
axis(axHm1_1,'off');
colormap(axHm1_1,colormapshow);
set(axHm1_1,'Clim',[-maxv maxv])
hold(axHm1_1,'on');
for i = 1:NROI
    plot([0.5,NROI+0.5],[0.5+(i-1),0.5+(i-1)],'color','k','linewidth',width1,'parent',axHm1_1);
    plot([0.5+(i-1),0.5+(i-1)],[0.5,NROI+0.5],'color','k','linewidth',width1,'parent',axHm1_1);
end
for i = 1:length(enhanind)
    plot([0.5,NROI+0.5],[0.5+enhanind(i),0.5+enhanind(i)],'color','k','linewidth',width2,'parent',axHm1_1);
    plot([0.5+enhanind(i),0.5+enhanind(i)],[0.5,NROI+0.5],'color','k','linewidth',width2,'parent',axHm1_1);
end

saveas(Hm1,[outdir,filesep,outname,'.fig'])
set(Hm1,'PaperPositionMode','manual');
set(Hm1,'PaperUnits','inch')
set(Hm1,'Paperposition',[1,1,POSm1(3)*3/300,POSm1(4)*3/300]);
print(Hm1,[outdir,filesep,outname,'.tif'],'-dtiff','-r300')
title(outname)
end