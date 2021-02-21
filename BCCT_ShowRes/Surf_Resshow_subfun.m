function Surf_Resshow_subfun(Parameter)
pathfiles = which('BCCT.m');
[pat nam ext] = fileparts(pathfiles);
Outdir = Parameter.Outputdir;

% background = [path,filesep,'templates',filesep];
FStype = {'fsaverage','fsaverage3','fsaverage4','fsaverage5','fsaverage6','fsaverage_sym'};

[surfp,ab] = SurfStatReadSurf({[pat,filesep,'templates',filesep,FStype{Parameter.FSTYPE},filesep,'lh.pial'],[pat,filesep,'templates',filesep,FStype{Parameter.FSTYPE},filesep,'rh.pial']});
[surfw,ab] = SurfStatReadSurf({[pat,filesep,'templates',filesep,FStype{Parameter.FSTYPE},filesep,'lh.white'],[pat,filesep,'templates',filesep,FStype{Parameter.FSTYPE},filesep,'rh.white']});
surfshow = surfp;
surfshow.coord = (surfp.coord+surfw.coord)/2;
save test2
Pvalthr = Parameter.PVAL;
if Parameter.PermTestlab==0 % not permutation test
    if Parameter.DATATYPE ==1
        ST_dat = SurfStatReadData(Parameter.LHRHfile);
        P_dat = SurfStatReadData(Parameter.LHRHfile_P);
        Pindex = P_dat<Pvalthr;
        Outshow{1} = ST_dat.*Pindex;
    elseif Parameter.DATATYPE==2
        ST_dat = SurfStatReadData(Parameter.LHRHfile);
        try 
            P_dat = SurfStatReadData(Parameter.LHRHfile_P);
            Pindex = P_dat>(1-Pvalthr);
            Outshow{1} = ST_dat.*Pindex;
        catch % only for NetFx2y
            Outshow{1} = ST_dat;
        end
        if Parameter.PermPLab
            Perm_dat = SurfStatReadData(Parameter.LHRHfile_Perm);
            Pindex = Perm_dat>(1-Pvalthr);
            Outshow{2} = ST_dat.*Pindex;
        end
    else
        ST_dat = SurfStatReadData(Parameter.LHRHfile);
        P_dat = SurfStatReadData(Parameter.LHRHfile_P);        
        Pindex = P_dat<Pvalthr;
        Outshow{1} = ST_dat.*Pindex;
        if Parameter.PermPLab
            Perm_dat = SurfStatReadData(Parameter.LHRHfile_Perm);
            Perm_dat2 = Perm_dat;
            Perm_dat2(Perm_dat>0.5) = 1-Perm_dat(Perm_dat>0.5);
            Pindex = Perm_dat2<Pvalthr;
            Outshow{2} = ST_dat.*Pindex;
        end
    end
else % permutation test
    ST_dat = SurfStatReadData(Parameter.PermLHRHfile);
    ind1 = find(ST_dat>(1-Pvalthr)&ST_dat<1);
    ind2 = find(ST_dat<Pvalthr&ST_dat>0);
    showmap = zeros(size(ST_dat));
    showmap(ind1) = ST_dat(ind1);
    showmap(ind2) = ST_dat(ind2)-1;
    Outshow{1} = showmap;
end


ColormapOut = AFNICOLORMAP(64);
% ColormapOutPos = ColormapOut(1:32,:);
% ColormapOutNeg = ColormapOut(33:64,:);
ColormapOut2 = [ColormapOut(1:32,:);0.8,0.8,0.8;ColormapOut(33:64,:)];

name = Parameter.name;
ex_name = {['P<',num2str(Pvalthr)],['PermP<',num2str(Pvalthr)]};
ex_name2 = {['Pless',num2str(Pvalthr)],['PermPless',num2str(Pvalthr)]};
for i = 1:length(Outshow)
    h=figure;
    SurfStatView(Outshow{i},surfshow,[name,ex_name{i}]);
    maxv = max(abs(Outshow{i}));
    if maxv<=0
        maxv = 0.5;
    end
    SurfStatColLim([-1,1]*maxv);
    SurfStatColormap(ColormapOut2);
    saveas(h,[Outdir,filesep,name,ex_name2{i},'.fig']);
    posnew = get(h,'pos');
    
%     set(h,'PaperPositionMode','manual');
%     set(h,'PaperUnits','inch')
%     set(h,'Paperposition',posnew);
    print(h,[Outdir,filesep,name,ex_name2{i},'.tif'],'-dtiff','-r300')
    
end
uiwait(msgbox('Finished'));

end