function [recon1_CS]=...
    XDGRASPRecon(recon1_nufft,param,lambdaT,lambdaS,idxGroup)

x = recon1_nufft;
scale = 1;
param.F = @(z)ForNUFFT_GPU(param,z)*scale;
param.Ft = @(z)AdjNUFFT_GPU(param,z,[size(x,1),size(x,2),size(x,3)])*scale;
scale = sqrt(1/(eps+abs(mean(vec(param.Ft(param.F(ones(size(x)))))))));
param.F = @(z)ForNUFFT_GPU(param,z)*scale;
param.Ft = @(z)AdjNUFFT_GPU(param,z,[size(x,1),size(x,2),size(x,3)])*scale;

recon1_CS = param.Ft(param.y);
normal = 1/max(abs(recon1_CS(:)));
param.y = param.y*normal;
recon1_CS = param.Ft(param.y);

param.TV_TempWeightT = max(abs(recon1_CS(:)))*lambdaT;
param.TV_TempWeightS = max(abs(recon1_CS(:)))*lambdaS;
param.TV_TempT = TV_Temp();
param.TV_TempS = TV_Temp_3D();
param.nite = 8;
param.display = 1;
tic

for n = 1 : 2
    n
    recon1_CS = CSL1NlCg_XDGrasp(recon1_CS,param);
    showDyImg(recon1_CS(:,:,:,[1,end]),idxGroup,300+n)
end
toc;
t_CS = toc;
t_CS = t_CS / 60