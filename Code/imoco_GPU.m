function [X] = imoco_GPU(param,dcf,traj,mr_img,smap,reg_field, TGV_lambda,ref,r,reg_field2)
% mc recon script
% Input:
%   fname: data file name base
% Output:
%   X:  imoco recon
% Xucheng Zhu, August 2019

% GPU version, Zekang Ding, 2020 07

addpath(genpath('../../pics'));

% motion resolved data shape[1,readout,spokes,sens,motion phase]
% data = readcfl_s([fname_base,'_datam']);
data = permute(param.y,[5,1,2,3,4]);
% data = single(data/max(abs(data(:))));
% motion resolved dcf shape[1,readout,spokes,1,motion phase]
% dcf = readcfl_s([fname_base,'_dcf2m']);
dcf = permute(dcf,[4,1,2,5,3]);
% motion resolved traj shape[3,readout,spokes,1,motion phase]
% traj =  readcfl_s([fname_base,'_trajm']);
traj = permute(traj,[1,2,3,5,4]);
% motion resolved recon shape[Xm,Ym,Zm,1,motion phase]
% mr_img = readcfl_s([fname_base,'_mrL']);
mr_img = permute(mr_img,[1,2,3,5,4]);
% motion resolved recon sensitivity[X,Y,Z,sens,1,1]
% smap = readcfl_s([fname_base,'_maps']);

IsizeL = size(mr_img);
m_ph = IsizeL(end);
s_datac = single(size(smap));

if nargin < 7
    TGV_lambda = .01;
end

nframe = m_ph;
coeff = 1:nframe;
coeff = (nframe - abs(coeff - ref)).^r;
weighting = nframe*coeff/sum(coeff);
weighting = sqrt(weighting);

for i = 1:nframe
    data(:,:,:,:,i) = data(:,:,:,:,i) * weighting(i);
end

%%%%%%%%%%

% width = 3;
% kb_1d = kaiser(2*width+1,13.9086);
% [kb_x,kb_y,kb_z] = meshgrid(kb_1d,kb_1d,kb_1d);
% kb_ksp = pad3d(kb_x.*kb_y.*kb_z,s_datac(1:3));
% kb_scale = (numel(kb_ksp)/((2*width)^3))*ifft3c(kb_ksp);
% 
% smap_kb =  smap./kb_scale;%wrong
smap_kb = smap;

% trick

% reg_field2 = inv_field(reg_field);

scale = 1;
E = GeneDygpu_NUFFTOperator(squeeze(traj),squeeze(dcf),smap_kb,8);
Afun = @(x)Forward(x,param,E,-reg_field2,weighting)*scale;
ATfun = @(y)Adjacent(squeeze(y),E,reg_field,weighting)*scale;
scale = sqrt(1/(eps+abs(mean(vec(ATfun(Afun(ones(size(smap(:,:,:,1))))))))));
Afun = @(x)Forward(x,param,E,-reg_field2,weighting)*scale;
ATfun = @(y)Adjacent(squeeze(y),E,reg_field,weighting)*scale;

% Afun = @(x)WGFSM1(x,traj,dcf,reg_field,smap_kb)*scale;
% ATfun = @(y)WGFSM1_H(y,traj,dcf,reg_field,smap_kb)*scale;
% scale = sqrt(1/(eps+abs(mean(vec(ATfun(Afun(ones(size(smap(:,:,:,1))))))))));
% Afun = @(x)WGFSM1(x,traj,dcf,reg_field,smap_kb)*scale;
% ATfun = @(y)WGFSM1_H(y,traj,dcf,reg_field,smap_kb)*scale;

Y = complex(double(zeros(size(data))));
Y_b = Y;
sigma = .25;
tau = .2;
rho = 5e-1;
datadcf = data;
% datadcf = data.*dcf;
X = ATfun(datadcf);
X_b = X;

% TGV parameters
params.lambda = TGV_lambda;
params.sigma = .25;
params.tau = .25;
params.alpha0 = .01;
params.alpha1 = .01;
params.nflag = 1;
TGV_prox = TGV(params);
tic;
theta = 1;
for iter_out = 1:20
    Y = (Y+sigma*(Afun(X_b)-datadcf))/(sigma + 1);
    X_1 = TGV_prox*(X-tau*ATfun(Y));
    X_b = X_1 + theta*(X_1-X);
    X = X_1;
    fprintf('Iter:%d, update_norm:%f.\n',iter_out,norm(X(:)-X_b(:))./norm(X(:)));
    if mod(iter_out,5)==0 || iter_out ==1
        showDyImg(X(:,:,:),[250 95 160],200+iter_out)
    end
end

% I_imoco = X;

mc_time = toc
% save([fname_base,'_imoco_pd',num2str(m_ph),'.mat'],'I_imoco','I_moco');
% I_sg = readcfl_s([fname_base,'_sg']);
% save([fname_base,'moco_pd.mat'],'X','mr_img','Ix','I_sg');%,'I_sg');
end

function iMfield = inv_field(Mfield)
iMfield = -Mfield;
for i = size(Mfield,4)
    iMfield(:,:,:,i,:) = -imwarp4(squeeze(Mfield(:,:,:,i,:)),-Mfield);
end
end

function z = imwarp4(x,y)
z = zeros(size(x));
for i = 1:size(y,5)
    z(:,:,:,i) = imwarp(x(:,:,:,i),y(:,:,:,:,i));
end
end


