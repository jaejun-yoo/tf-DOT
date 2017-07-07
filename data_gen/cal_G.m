function G_r1_r2 = cal_G(k,D0,r1X,r1Y,r1Z,r2X,r2Y,r2Z)
r1_coordi = [r1Y(:),r1X(:), r1Z(:)];
r2_coordi = [r2Y(:),r2X(:), r2Z(:)];
dMat = pdist2((r1_coordi),(r2_coordi),'euclidean'); % |r_1 - r_2|
% G_r1_r2 = bsxfun(@rdivide,exp(-1*k*dMat),(4*pi*D0*dMat)); 
G_r1_r2 = exp(-1*k*dMat)./(4*pi*D0*dMat);
end