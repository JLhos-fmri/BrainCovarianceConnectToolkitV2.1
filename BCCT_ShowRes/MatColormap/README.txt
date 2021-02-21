
load FCmat % load Functional conn matrix
load mycolormap % load colormap from brewermap.m ( 2015 Stephen Cobeldick)


figure(1),imagesc(ROICorrelation)
colormap(flipud(RdYlBu_map))  % RdYlBu


figure(2),imagesc(ROICorrelation)
colormap(flipud(RdBu_map))   % RdBu


figure(3),imagesc(ROICorrelation)
colormap(flipud(Spectral_map))   % Spectral


figure(4),imagesc(ROICorrelation)
colormap(esa(256))


figure(5),imagesc(ROICorrelation)
colormap(redblue(256))


% figure(4),imagesc(ROICorrelation)
%  colormap(darkb2r(-1,1))