function ksp = ForNUFFT_GPU(param,img)

%DingZekang 20200711 at UIH
% kdata_cs  RO LPE nc nt

FT = param.E;
[RO1,LPE,LCH_PCA,tt] = size(param.y);

ksp = zeros(RO1*LPE,LCH_PCA,tt);

% tmppath = pwd;
for i = 1:tt
    single_img = img(:,:,:,i);
    ksp(:,:,i) = FT{i}*single_img;
end
ksp = reshape(ksp,[RO1,LPE,LCH_PCA,tt]);