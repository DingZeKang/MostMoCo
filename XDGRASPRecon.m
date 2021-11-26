function [recon1_CS]=...
    XDGRASPRecon(recon1_nufft,param,lambdaT,lambdaS,idxGroup)

x = recon1_nufft;
scale = 1;
param.F = @(z)ForNUFFT_GPU(param,z)*scale;
param.Ft = @(z)AdjNUFFT_GPU(param,z,[size(x,1),size(x,2),size(x,3)])*scale;
scale = sqrt(1/(eps+abs(mean(vec(param.Ft(param.F(ones(size(x)))))))));
param.F = @(z)ForNUFFT_GPU(param,z)*scale;
param.Ft = @(z)AdjNUFFT_GPU(param,z,[size(x,1),size(x,2),size(x,3)])*scale;

%param.y%0.0155
recon1_CS = param.Ft(param.y);%0.0063
normal = 1/max(abs(recon1_CS(:)));%157
param.y = param.y*normal;%2.4499
recon1_CS = param.Ft(param.y);%1.0

param.TV_TempWeightT = max(abs(recon1_CS(:)))*lambdaT;%for 16 times undersampling,1,0.5 is wrong,0.1 ok,0.1 better than 0.02
param.TV_TempWeightS = max(abs(recon1_CS(:)))*lambdaS;
param.TV_TempT = TV_Temp();
param.TV_TempS = TV_Temp_3D();
param.nite = 8;
param.display = 1;
tic

for n = 1 : 2% 
    n
    recon1_CS = CSL1NlCg_XDGrasp(recon1_CS,param);
    showDyImg(recon1_CS(:,:,:,[1,end]),idxGroup,300+n)
end
toc;
t_CS = toc;
t_CS = t_CS / 60