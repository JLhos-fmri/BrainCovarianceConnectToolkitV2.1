function ColormapOut = AFNICOLORMAP(k)
if mod(k,2)
    k = k+1;
end
mk = k/2;
ColorMapOrig = [1,1,0;1,0.8,0;1,0.6,0;1,0.4118,0;1,0.2667,0;1,0,0;0,0,1;0,0.2667,1;0,0.4118,1;0,0.6,1;0,0.8,1;0,1,1;];
ColorMapMid = ColorMapOrig(:,2);
switch k
    case 2,
        ColorMap=[0,1,1;
            1,0.8,0;];
    case 4,
        ColorMap=[0,1,1;
            0,0.4118,1;
            1,0.2667,0;
            1,0.8,0;];
    case 6,
        ColorMap=[0,1,1;
            0,0.6,1;
            0,0.2667,1;
            1,0,0;
            1,0.4118,0;
            1,0.8,0;];
    case 8,
        ColorMap=[0,1,1;
            0,0.8,1;
            0,0.4118,1;
            0,0.2667,1;
            1,0.2667,0;
            1,0.4118,0;
            1,0.6,0;
            1,0.8,0;];
    case 10,
        ColorMap=[0,1,1;
            0,0.8,1;
            0,0.6,1;
            0,0.4118,1;
            0,0.2667,1;
            1,0,0;
            1,0.2667,0;
            1,0.4118,0;
            1,0.6,0;
            1,0.8,0;];
    case 12
        ColorMap = ColorMapOrig(end:-1:1,:);
    otherwise
        ColorMapMid2 = interp1(1:12,ColorMapMid,1:11/(k-1):12);
        ColorMap = [[zeros(mk,1);ones(mk,1)],ColorMapMid2',[ones(mk,1);zeros(mk,1)]];
end
ColormapOut = ColorMap;