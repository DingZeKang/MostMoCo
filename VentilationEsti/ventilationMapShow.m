close all
clc
clear
path = '/home/dzk/DingZekang/DataAnalysis/20210508New/';
name = {'DZK','HZ','LB','LQK','LRH','QSH','SHJ','WRK','YK','ZZX'};
subject = 10;
siz = [496,348,372];
mask = zeros([siz,subject]);
Vmap1 = zeros([siz,subject]);
Vmap2 = zeros([siz,subject]);
Vmap3 = zeros([siz,subject]);
Vmap4 = zeros([siz,subject]);
Fmap1 = zeros([siz,5,subject]);
Fmap2 = zeros([siz,5,subject]);
Fmap3 = zeros([siz,5,subject]);
Fmap4 = zeros([siz,5,subject]);
AligImg_XD = zeros([siz,5,subject]);
AligImg_iMoCo = zeros([siz,5,subject]);
AligImg_Most = zeros([siz,5,subject]);
AligImg_SENSE = zeros([siz,5,subject]);
for i = 1:subject
    load([path,name{i},'/mask3D2.mat'])
    mask(:,:,:,i) = BW;
    load([path,name{i},'/Rep8AligImg_XD.mat']);
    AligImg_XD(:,:,:,:,i) = AligImg;
    load([path,name{i},'/Rep8AligImg_iMoCo.mat']);
    AligImg_iMoCo(:,:,:,:,i) = AligImg;
    load([path,name{i},'/Rep8AligImg_Most.mat']);
    AligImg_Most(:,:,:,:,i) = AligImg;
end
mask = logical(mask);

sigma = 3;
range = 2:4;
for i = 1:subject
    i
    idx = find(mask(:,:,:,i) ~= 0);
    zidx = ceil(idx/(496*348));
    yidx = ceil((idx-(zidx-1)*496*348)/496);
    xidx = idx-(zidx-1)*496*348-(yidx-1)*496;
    Fmap1(:,:,:,:,i) = calven_FT(AligImg_XD(:,:,:,:,i),sigma,xidx,yidx,zidx);
    Fmap2(:,:,:,:,i) = calven_FT(AligImg_iMoCo(:,:,:,:,i),sigma,xidx,yidx,zidx);
    Fmap3(:,:,:,:,i) = calven_FT(AligImg_Most(:,:,:,:,i),sigma,xidx,yidx,zidx);
%     Fmap4(:,:,:,:,i) = calven_FT(AligImg_SENSE(:,:,:,:,i),sigma,xidx,yidx,zidx);
end
Fmap1 = abs(Fmap1);Fmap2 = abs(Fmap2);Fmap3 = abs(Fmap3);Fmap4 = abs(Fmap4);
save(['maskFrame3All.mat'],'mask','-v7.3');
save(['Rep8Fmap1_sig',num2str(sigma),'_fit6.mat'],'Fmap1','-v7.3');
save(['Rep8Fmap2_sig',num2str(sigma),'_fit6.mat'],'Fmap2','-v7.3');
save(['Rep8Fmap3_sig',num2str(sigma),'_fit6.mat'],'Fmap3','-v7.3');

%%% another process
load('Rep8Fmap1_sig3_fit6.mat');load('Rep8Fmap2_sig3_fit6.mat');load('Rep8Fmap3_sig3_fit6.mat');load('Fmap4_sig2_fit6.mat');
load('maskFrame3All.mat');
masker = mask;
range = 2:4;

for i = 1:subject
    ii = i
    Vmap1(:,:,:,i) = calven(Fmap1(:,:,:,:,i),range);
    Vmap2(:,:,:,i) = calven(Fmap2(:,:,:,:,i),range);
    Vmap3(:,:,:,i) = calven(Fmap3(:,:,:,:,i),range);
    Vmap4(:,:,:,i) = calven(Fmap4(:,:,:,:,i),range);
end

traidx = 195;
coridx = 105;
sagidx = 185;
idxGroup = [traidx,coridx,sagidx];

i = 1;
tmp = cat(5,Vmap4,Vmap1,Vmap2,Vmap3);
showDyImgMF(squeeze(tmp(:,:,:,10,:)),idxGroup,50)

tmpvect0 = Vmap4(masker);
tmpvect0 = tmpvect0(~isnan(tmpvect0));
meanVen0 = mean(tmpvect0);
stdVen0 = std(tmpvect0);

tmpvect1 = Vmap1(masker);
tmpvect1 = tmpvect1(~isnan(tmpvect1));
meanVen1 = mean(tmpvect1);
stdVen1 = std(tmpvect1);

tmpvect2 = Vmap2(masker);
tmpvect2 = tmpvect2(~isnan(tmpvect2));
meanVen2 = mean(tmpvect2);
stdVen2 = std(tmpvect2);

tmpvect3 = Vmap3(masker);
tmpvect3 = tmpvect3(~isnan(tmpvect3));
meanVen3 = mean(tmpvect3);
stdVen3 = std(tmpvect3);

[h1,p1]=ttest(tmpvect1,tmpvect3,0.1)
[h2,p2]=ttest(tmpvect2,tmpvect3,0.1)

a = 1;

close all
Vmap_all = cat(5,Vmap4,Vmap1,Vmap2,Vmap3);
for i = 1:10
    tmpimgAll = squeeze(Vmap_all(:,100,:,i,1:end));
    tmpimg = permute(tmpimgAll(:,:,:),[2,1,3]);
    tmpimg = reshape(tmpimg,size(tmpimg,1),size(tmpimg,2)*size(tmpimg,3));
    figure,imshow((tmpimg),[0,0.3],'border','tight'),colormap Jet,%colorbar;
end