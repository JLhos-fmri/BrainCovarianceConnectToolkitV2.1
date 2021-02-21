function dynamicBC_Reslice(PI,PO,NewVoxSize,hld,TargetSpace)
%   PI - input filename
%   PO - output filename
%   NewVoxSize - 1x3 matrix of new vox size.
%   hld - interpolation method. 0: Nearest Neighbour. 1: Trilinear.
%   TargetSpace - Define the target space. 'ImageItself': defined by the image itself (corresponds  to the new voxel size); 'XXX.img': defined by a target image 'XXX.img' (the NewVoxSize parameter will be discarded in such a case).
if nargin<=4
    TargetSpace='ImageItself';
end

if ~strcmpi(TargetSpace,'ImageItself')
    headIN = spm_vol(TargetSpace) ;
    headIN = headIN(1);
    dataIN = spm_read_vols(headIN);
    mat=headIN.mat;
    dim=headIN.dim;
else
    headIN = spm_vol(PI) ;
    headIN = headIN(1);
    dataIN = spm_read_vols(headIN);
    origin=headIN.mat(1:3,4);
    origin=origin+[headIN.mat(1,1);headIN.mat(2,2);headIN.mat(3,3)]-[NewVoxSize(1)*sign(headIN.mat(1,1));NewVoxSize(2)*sign(headIN.mat(2,2));NewVoxSize(3)*sign(headIN.mat(3,3))];
    origin=round(origin./NewVoxSize').*NewVoxSize';
    mat = [NewVoxSize(1)*sign(headIN.mat(1,1))                 0                                   0                       origin(1)
        0                         NewVoxSize(2)*sign(headIN.mat(2,2))              0                       origin(2)
        0                                      0                      NewVoxSize(3)*sign(headIN.mat(3,3))  origin(3)
        0                                      0                                   0                          1      ];

    dim=(headIN.dim-1).*diag(headIN.mat(1:3,1:3))';
    dim=floor(abs(dim./NewVoxSize))+1;
end
VI          = spm_vol(PI);
VI   = VI(1);
VO          = VI;
VO.fname    = deblank(PO);
VO.mat      = mat;
VO.dim(1:3) = dim;

[x1,x2] = ndgrid(1:dim(1),1:dim(2));
d     = [hld*[1 1 1]' [1 1 0]'];
C = spm_bsplinc(VI, d);
v = zeros(dim);
for x3 = 1:dim(3),
    [tmp,y1,y2,y3] = getmask(inv(mat\VI.mat),x1,x2,x3,VI.dim(1:3),[1 1 0]');
    v(:,:,x3)      = spm_bsplins(C, y1,y2,y3, d);
end;
VO = spm_write_vol(VO,v);
end

function [Mask,y1,y2,y3] = getmask(M,x1,x2,x3,dim,wrp)
tiny = 5e-2; % From spm_vol_utils.c
y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
Mask = logical(ones(size(y1)));
if ~wrp(1), Mask = Mask & (y1 >= (1-tiny) & y1 <= (dim(1)+tiny)); end;
if ~wrp(2), Mask = Mask & (y2 >= (1-tiny) & y2 <= (dim(2)+tiny)); end;
if ~wrp(3), Mask = Mask & (y3 >= (1-tiny) & y3 <= (dim(3)+tiny)); end;
return;
end