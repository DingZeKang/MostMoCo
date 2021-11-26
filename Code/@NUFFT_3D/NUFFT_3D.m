function  res = NUFFT_3D(k,w,b1)

% NUFFT_3D 
% Shuo Li, SJTU, 2019
% based on:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multicoil NUFFT operator
% Based on the NUFFT toolbox from Jeff Fessler and the single-coil NUFFT
% operator from Miki Lustig
% Input
% k: k-space trajectory        /3,samples,total_spokes,frames/
% w: density compensation      /samples,total_spokes,frames/
% b1: coil sensitivity maps    /res,res,slices/  for one-coil
%
% Li Feng & Ricardo Otazo, NYU, 2012
    
Nd = size(b1(:,:,:,1));%Nd=[256 256 256]
Jd = [3,3,3];
Kd = floor([Nd*1.5]);%kd=[384 384 384]
n_shift = Nd/2;%n_shift=[128 128 128]
for tt = 1:size(k,4)
	kk = k(:,:,:,tt);%k=tu
	om = [kk(1,:)',kk(2,:)',kk(3,:)']*2*pi;
	res.st{tt} = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
end
res.adjoint = 0;
res.imSize = size(b1(:,:,:,1));
res.imSize2 = 2*size(b1(:,:,:,1));
res.dataSize = [size(k,2),size(k,3)];
res.w = sqrt(w);
res.b1 = b1;
res = class(res,'NUFFT_3D');

