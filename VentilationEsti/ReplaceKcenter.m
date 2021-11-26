% kspace center replacement
% clear
clc
close all
load('recon_nufft.mat')
recon1_nufft = recon1_nufft/max(abs(recon1_nufft(:)));
load('XDGRASPImage.mat');
recon_XDGRASP = recon_XDGRASP/max(abs(recon_XDGRASP(:)));
load('iMoCoImage.mat');
DysImage = DysImage/max(abs(DysImage(:)));
load('MostMoCoImage.mat');%
DyImage = DyImage/max(abs(DyImage(:)));

traidx = 250;%porcine
coridx = 95;%porcine
sagidx = 160;%porcine
idxGroup = [traidx,coridx,sagidx];
idxGroup = [120,140,150];%human

siz = size(DyImage);
nframe = siz(4);

showDyImg(recon1_nufft(:,:,:,[1,nframe]),round(idxGroup),601)
showDyImg(recon_XDGRASP(:,:,:,[1,nframe]),round(idxGroup),602)
showDyImg(DysImage(:,:,:,[1,nframe]),round(idxGroup),603)
showDyImg(DyImage(:,:,:,[1,nframe]),round(idxGroup),604)

cx = round(siz(1)/2); cy = round(siz(2)/2); cz = round(siz(3)/2);

originKsp = fft4(recon1_nufft);
showDyImg(originKsp(:,:,:,[1,nframe]),round([cx,cy,cz]),600)

Reprecon_XDGRASP = replaceFunc(recon_XDGRASP,originKsp,E,Et);
RepDysImage = replaceFunc(DysImage,originKsp,E,Et);
RepDyImage = replaceFunc(DyImage,originKsp,E,Et);

showDyImg(Reprecon_XDGRASP(:,:,:,[1,nframe]),round(idxGroup),605)
showDyImg(RepDysImage(:,:,:,[1,nframe]),round(idxGroup),606)
showDyImg(RepDyImage(:,:,:,[1,nframe]),round(idxGroup),607)

max(abs(Reprecon_XDGRASP(:)))
max(abs(RepDysImage(:)))
max(abs(RepDyImage(:)))

Reverse = 1;% from End Expiratio to End Inspiration
if Reverse
%      recon1_nufft = recon1_nufft(:,:,:,end:-1:1);
     Reprecon_XDGRASP = Reprecon_XDGRASP(:,:,:,end:-1:1);
     RepDysImage = RepDysImage(:,:,:,end:-1:1);
     RepDyImage = RepDyImage(:,:,:,end:-1:1);
end

M2 = Reprecon_XDGRASP; M3 = RepDysImage; M4 = RepDyImage;
showDyImg(M2(:,:,:,[1,nframe]),round([idxGroup]),500)

