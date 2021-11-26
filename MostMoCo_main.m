%Motion State weighted Motion Compensation Dynamic Pulmonary MR

% ZekangDing, Shanghai Jiao Tong University
%20211125

% main program
clc
clear 
close all
addpath(genpath('Code'));
addpath(genpath('dataExample'));
setenv('CUDA_VISIBLE_DEVICES','0');
cardIdx = 1;
gpuDevice(cardIdx);

%%%%% before recon, edit the recon parameters!!!!%%%%%%%%%%%
ratio = 1;
pathname1='reconResult/';%data saving path

UsingGroupReg = 0;
ThrowSeg = 0;

nframe = 5;%5          

traidx = 110;
coridx = 105;
sagidx = 120;
idxGroup = [traidx,coridx,sagidx];

% Performing following steps to run this program:
% 1. Read Rawdata, SI navigator data(kdata_SI), and sequence parameter(seqParam);
% 2. Calculate trajectory(Crds) and density compensation(DCF);
% 3. Apply coil compression on rawdata, get kdata_cc;

load('kdata_cc.mat');load('Crds.mat');load('DCF.mat');
load('kdata_SI.mat');load('seqParam.mat');
UsingSINav = 1;

numPoint = round(ratio*size(kdata_cc,1));

[imgNufft,imgLNufft,sens,sensL,k,w]=bartReconSens(kdata_cc,DCF,Crds,idxGroup,ratio);

img_size = size(imgNufft);
lowimg_size = size(imgLNufft);
Imgratio = lowimg_size(1)/img_size(1);

if UsingSINav
    TrueIndex = RespirationEsti(kdata_SI,ThrowSeg,seqParam);
else
    kc0 = Extractk0(FIDrawdata,seqParam);
    TrueIndex = RespirationEsti_KC0((kc0),ThrowSeg,seqParam);
end

kdata_bart = kdata_cc/max(abs(kdata_cc(:)));
kdata_bart=kdata_bart(:,TrueIndex,:);
k_bart=k(:,:,TrueIndex);
w_bart=w(:,TrueIndex);
[k_bart,w_bart,kdata_cs] = data_sorting(k_bart,w_bart,kdata_bart,nframe);
kdata_cs = permute(kdata_cs,[7,1,2,3,6,4,5]);
k_bart = permute(k_bart,[1,2,3,7,6,4,5]);
w_bart = permute(w_bart,[7,1,2,6,5,3,4]);
w_bart = w_bart/max(abs(w_bart(:)));
w_cs2 = sqrt(w_bart);

writecfl('k_cs',k_bart);
writecfl('w_cs',w_bart);
writecfl('w_cs2',w_cs2);
writecfl('k_csL',k_bart(:,1:numPoint,:,:,:,:,:));
writecfl('w_cs2L',w_cs2(:,1:numPoint,:,:,:,:,:));

kdatadcf = bart('fmac',kdata_cs,w_bart);
img = bart('nufft -a',k_bart(:,1:numPoint,:,:,:,:,:),kdatadcf(:,1:numPoint,:,:,:,:,:));
lowres_nufft = squeeze(bart('fmac -C -s 8',img,sensL));
showDyImg(squeeze(lowres_nufft(:,:,:,[1,nframe])),round(idxGroup*Imgratio),700)

lowres_CS = squeeze(bart('pics -C 20 -i 20 -R T:7:0:0.005 -p w_cs2L -t k_csL',kdata_cs(:,1:numPoint,:,:,:,:,:),sensL));
showDyImg(squeeze(lowres_CS(:,:,:,[1,nframe])),round(idxGroup*Imgratio),701)

[recon_nufft,param,k_cs,w_cs] = Dynamic_NUFFT(squeeze(kdata_cs.*w_cs2),squeeze(k_bart/max(k_bart(:))/2),squeeze(w_bart),img_size,sens,nframe,8);
showDyImg(recon_nufft(:,:,:,[1,nframe]),round(idxGroup),604)
clear kdata_bart FIDrawdata img kdatadcf

save([pathname1,'recon_nufft.mat'],'recon_nufft','-v7.3');

% XD-UTE recon
[recon_XDGRASP]=XDGRASPRecon(recon_nufft,param,0.05,0.01,idxGroup);
save([pathname1,'XDGRASPImage.mat'],'recon_XDGRASP','-v7.3');

% MostMoCo recon
ref = round(nframe/2);

%Image Registration
Img_L = abs(lowres_CS)./max(abs(lowres_CS(:)))+eps;
nbin = size(Img_L,4);
sizeH = size(recon_nufft(:,:,:,1));
sizeL = size(Img_L(:,:,:,1));
downsampling = sizeH./sizeL;

[B_L,F_L] = registrationPandG(ref,nbin,Img_L,sizeL,UsingGroupReg,seqParam);
gpuDevice(cardIdx);


B = zeros([sizeH,3,nbin,nbin]);
F = zeros([sizeH,3,nbin]);
[X1,Y1,Z1] = meshgrid(linspace(1,sizeL(2),sizeL(2)),...
    linspace(1,sizeL(1),sizeL(1)),...
    linspace(1,sizeL(3),sizeL(3)));
[Xq,Yq,Zq] = meshgrid(linspace(1,sizeL(2),sizeH(2)),...
    linspace(1,sizeL(1),sizeH(1)),...
    linspace(1,sizeL(3),sizeH(3)));
for ref = 1:nbin
for i = 1:nbin
    for j = 1:3       
        B(:,:,:,j,i,ref) = downsampling(j)*interp3(X1,Y1,Z1,B_L(:,:,:,j,i,ref),Xq,Yq,Zq,'cubic');
        F(:,:,:,j,i) = downsampling(j)*interp3(X1,Y1,Z1,F_L(:,:,:,j,i),Xq,Yq,Zq,'cubic');
    end
end
end
showDyImgMF(B(:,:,:,3,1,3),idxGroup,200)
clear X1 Y1 Z1 Xq Yq Zq

x = recon_nufft;
scale = 1;
E = @(z)ForNUFFT_GPU(param,z)*scale;
Et = @(z)AdjNUFFT_GPU(param,z,[size(x,1),size(x,2),size(x,3)])*scale;
scale = sqrt(1/(eps+abs(mean(vec(Et(E(ones(size(x)))))))));
E = @(z)ForNUFFT_GPU(param,z)*scale;
Et = @(z)AdjNUFFT_GPU(param,z,[size(x,1),size(x,2),size(x,3)])*scale;

recon1_CS = Et(param.y);
normal = 1/max(abs(recon1_CS(:)));
param.y = param.y*normal;
clear recon1_CS

mTVs = TVs(3);
TVt = TV();
sizeI2 = size(recon_nufft);

rho = 1;
rho1 = 1;

vec = @(z)z(:);
A = E;At = Et;
d = param.y;
Trans = Registra(nbin,B,F);
weightingCoef = 8;
mTVx = TVma(nbin,B,F,weightingCoef);

OutIter = 3;
Iter = 0;

% init
sizeL = sizeH;
downsampling = sizeH./sizeL;

x_k0 = At(d);
z_k0 = zeros(size(mTVx*x_k0));
u_k0 = zeros(size(z_k0));
z1_k0 = zeros(size(mTVs*x_k0));
u1_k0 = zeros(size(z1_k0));

lambda1 = 0.05*max(abs(x_k0(:)));
lambda2 = 0.01*max(abs(x_k0(:)));

while(Iter<OutIter)
    
    mTVx = TVma(nbin,B,F,weightingCoef);
    Trans = Registra(nbin,B,F);
    
    iter = 0;
    tic;
    while(iter<5)
        Iter,iter
        iter = iter+1;
        
        % temporal L1 normalization
        if rho
            if 1 %TVt of T_theta(x)
          
                z_k1 = wthresh(mTVx*(x_k0)+u_k0/rho,'s',lambda1/rho);
               
            else %TVt of operation(x)
                
            end
            test = squeeze(reshape(z_k1,sizeI2));
            showDyImg(test(:,:,:,[1,nframe]),idxGroup,900)
        end

        % spatial L1 normalization
        if rho1 %spatial regularization
            if 1 %spatial TV
                
                z1_k1 = wthresh(mTVs*(x_k0)+u1_k0/rho1,'s',lambda2/rho1);
          
            else %spatial wavelet
                
                %To do
                
            end
        end

        % L2 opt conjugate gradient descent
        cg_iterM = 30;
        tol = 1e-3;
        x_k1 = conj_grad_x_MotionAverage(A,At,x_k0,d,rho,mTVx,z_k1,u_k0,rho1,mTVs,z1_k1,u1_k0,tol,cg_iterM);
        
        test = squeeze(reshape(x_k1,sizeI2));
        showDyImg(test(:,:,:,[1,nframe]),idxGroup,902)
        
        % dual update
        if rho, u_k1 = u_k0 + rho*(mTVx*(x_k1) -z_k1); end
        if rho1, u1_k1 = u1_k0 + rho1*(mTVs*(x_k1) -z1_k1); end

        converg = norm(x_k1(:)-x_k0(:))/norm(x_k0(:))
        if converg <= 1e-5, break; end
        
        % all update
        z_k0 = z_k1; 
        z1_k0 = z1_k1;
        x_k0 = x_k1;
        u_k0 = u_k1;
        u1_k0 = u1_k1;
        rho = rho*1.0;
        rho1 = rho1*1.0;
        
    end
    % L1 wavelet
    % data consistancy
    toc
    
    tmpImg = squeeze(reshape(x_k0,sizeI2));
    showDyImg(tmpImg(:,:,:,[1,nframe]),idxGroup,1000+Iter)
    
    Img = x_k0;
    Img = squeeze(abs(Img)./max(abs(Img(:))))+eps;
    
    if Iter<OutIter-1      
        
        [X1,Y1,Z1] = meshgrid(linspace(1,sizeH(2),sizeH(2)),...
            linspace(1,sizeH(1),sizeH(1)),...
            linspace(1,sizeH(3),sizeH(3)));
        [Xq,Yq,Zq] = meshgrid(linspace(1,sizeH(2),sizeL(2)),...
            linspace(1,sizeH(1),sizeL(1)),...
            linspace(1,sizeH(3),sizeL(3)));
        for i = 1:nbin                   
            Img_tmp(:,:,:,i) = interp3(X1,Y1,Z1,Img(:,:,:,i),Xq,Yq,Zq,'cubic');            
        end
        [B_L,F_L] = registrationPandG(ref,nbin,Img_tmp,sizeL,UsingGroupReg,seqParam);
        gpuDevice(cardIdx);
        
        [X1,Y1,Z1] = meshgrid(linspace(1,sizeL(2),sizeL(2)),...
            linspace(1,sizeL(1),sizeL(1)),...
            linspace(1,sizeL(3),sizeL(3)));
        [Xq,Yq,Zq] = meshgrid(linspace(1,sizeL(2),sizeH(2)),...
            linspace(1,sizeL(1),sizeH(1)),...
            linspace(1,sizeL(3),sizeH(3)));
        for ref = 1:nbin
        for i = 1:nbin
            for j = 1:3       
                B(:,:,:,j,i,ref) = downsampling(j)*interp3(X1,Y1,Z1,B_L(:,:,:,j,i,ref),Xq,Yq,Zq,'cubic');
                F(:,:,:,j,i) = downsampling(j)*interp3(X1,Y1,Z1,F_L(:,:,:,j,i),Xq,Yq,Zq,'cubic');
            end
        end
        end
        
        showDyImgMF(B(:,:,:,3,1,3),idxGroup,200)
    end
    clear X1 Y1 Z1 Xq Yq Zq
    
    Iter = Iter + 1;    
end
DyImage = squeeze(x_k0);
showDyImg(DyImage(:,:,:,[1,nframe]),idxGroup,400)
save([pathname1,'MostMoCoImage.mat'],'DyImage','-v7.3');

%iMoCo
ratio = 0.75;
numPoint = round(ratio*size(kdata_cc,1));

[imgNufft,imgLNufft,sens,sensL,k,w]=bartReconSens(kdata_cc,DCF,Crds,idxGroup,ratio);
img_size = size(imgNufft);
lowimg_size = size(imgLNufft);
Imgratio = lowimg_size(1)/img_size(1);
writecfl('k_csL',k_bart(:,1:numPoint,:,:,:,:,:));
writecfl('w_cs2L',w_cs2(:,1:numPoint,:,:,:,:,:));

kdatadcf = bart('fmac',kdata_cs,w_bart);
img = bart('nufft -a',k_bart(:,1:numPoint,:,:,:,:,:),kdatadcf(:,1:numPoint,:,:,:,:,:));
lowres_nufft = squeeze(bart('fmac -C -s 8',img,sensL));
showDyImg(squeeze(lowres_nufft(:,:,:,[1,nframe])),round(idxGroup*Imgratio),700)

lowrecon1_CS = squeeze(bart('pics -C 20 -i 20 -R T:7:0:0.005 -p w_cs2L -t k_csL',kdata_cs(:,1:numPoint,:,:,:,:,:),sensL));%0.0001
showDyImg(squeeze(lowrecon1_CS(:,:,:,[1,nframe])),round(idxGroup*Imgratio),701)

mr_img = squeeze(lowrecon1_CS);
mr_img = abs(squeeze(mr_img)./max(abs(mr_img(:))));
IsizeL = size(mr_img);
m_ph = IsizeL(end);
IsizeL = IsizeL(1:3);
Isize = size(sens);
Isize = Isize(1:3);
mscale = Isize./IsizeL;
mask = ones(IsizeL);
% estimate motion state
[reg_field_L,~] = registrationPandG(ref,nbin,mr_img,IsizeL,UsingGroupReg,seqParam);
gpuDevice(cardIdx);
DysImage = zeros([Isize,nframe]);
for ref = 1:nframe
    ref
    reg_field = zeros([Isize(1:3),3,m_ph]);
    reg_field2 = zeros([Isize(1:3),3,m_ph]);
    % reg_field update  
   
    for i = 1:m_ph 
        for j = 1:3
            reg_field(:,:,:,j,i) = imresize3(reg_field_L(:,:,:,j,i,ref).*mask*mscale(j),Isize);
            reg_field2(:,:,:,j,i) = imresize3(reg_field_L(:,:,:,j,ref,i).*mask*mscale(j),Isize);
        end
    end

    sImage = imoco_GPU(param,(w_cs),k_cs,squeeze(lowrecon1_CS),sens,reg_field,0.1,1,0,reg_field2);

    showDyImg(sImage,idxGroup,611)
    DysImage(:,:,:,ref) = sImage;
    
end
showDyImg(DysImage(:,:,:,[1,nframe]),idxGroup,612)
filename=['iMoCoImage.mat'];
save([pathname1,filename],'DysImage','-v7.3');

