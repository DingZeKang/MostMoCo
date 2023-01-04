function img = AdjNUFFT_GPU(param,kdata_cs,img_size)

%DingZekang 20200711 at UIH
% kdata_cs  RO LPE nc nt

FT = param.E;
[RO1,LPE,LCH_PCA,tt] = size(kdata_cs);

kdata_cs = reshape(kdata_cs,[RO1*LPE,LCH_PCA,tt]);

img = zeros([img_size,tt]);

for i = 1:tt
    kdata = kdata_cs(:,:,i);
    img(:,:,:,i) = FT{i}'*(kdata);
end