function MSSIM = SSIMComp(x,y)

x = abs(x(:))/max(abs(x(:)))*255;
y = abs(y(:))/max(abs(y(:)))*255;

ux = mean(x); uy = mean(y);
sigmx = std(x); sigmy = std(y);
tmp = cov(x,y);
sigmxy = tmp(1,2);

k1 = 0.01; k2 = 0.03;
c1 = (k1*255)^2; c2 = (k2*255)^2;

MSSIM = (2*ux*uy+c1)*(2*sigmxy+c2)/(ux^2+uy^2+c1)/(sigmx^2+sigmy^2+c2);