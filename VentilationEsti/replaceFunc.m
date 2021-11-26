function RepImg = replaceFunc(Img,originKsp,E,Et)

siz = size(Img);
cx = round(siz(1)/2); cy = round(siz(2)/2); cz = round(siz(3)/2);
x = 8; y = 8; z = 8;

reconKsp = fft4(Img);
ReplaceKsp = reconKsp;
ReplaceKsp(cx-x/2+1:cx+x/2,cy-y/2+1:cy+y/2,cz-z/2+1:cz+z/2,:) = originKsp(cx-x/2+1:cx+x/2,cy-y/2+1:cy+y/2,cz-z/2+1:cz+z/2,:);
RepImg = ifft4(ReplaceKsp);

