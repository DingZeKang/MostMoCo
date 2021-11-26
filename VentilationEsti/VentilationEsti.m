
setenv('CUDA_VISIBLE_DEVICES','1');
cardIdx = 1;
gpuDevice(cardIdx);
traidx = 195;
coridx = 105;
sagidx = 185;
idxGroup = [traidx,coridx,sagidx];
siz = size(DyImage);
frame = siz(4);
load('mask3D2.mat');
idx = find(BW ~= 0);
zidx = floor(idx/(496*348))+1;
yidx = floor((idx-(zidx-1)*496*348)/496)+1;
xidx = idx-(zidx-1)*496*348-(yidx-1)*496;
ref = 3;
for numMet = 1:3
    if numMet == 1
        sMethod = 'XD';DyImage = M2;numMet
    elseif numMet == 2
        sMethod = 'iMoCo';DyImage = M3;numMet
    else
        sMethod = 'Most';DyImage = M4;numMet
    end
%     sMethod = 'Most';%iMoCo,XD,Most
%     DyImage = M4;
    B = zeros([siz(1:3),3,siz(4)]);
    AligImg = zeros(siz);
    Img = abs(DyImage/max(abs(DyImage(:))));
    Img = abs(imgauss4d(Img,.5));
    fixedGPU = gpuArray(Img(:,:,:,ref));
    fixedHist = imhist(fixedGPU);
    for i = [1:5]
        movingGPU = gpuArray(Img(:,:,:,i));
        movingGPU = histeq(movingGPU,fixedHist);

        [B_GPU, movingRegGPU] = imregdemons(movingGPU,fixedGPU,[200,100,50,25]*4,'PyramidLevels',4,'DisplayWaitbar',false,'AccumulatedFieldSmoothing',1.0);  %1.5 better
        B(:,:,:,:,i) = gather(B_GPU);
%         AligImg(:,:,:,i) = gather(movingRegGPU);%wrong, 'histeq' operation
%         changed the original intensity of moving image;
        AligImg(:,:,:,i) = imwarp(Img(:,:,:,i),B(:,:,:,:,i));
    end
    gpuDevice(cardIdx);
    showDyImg(AligImg(:,:,:,[1,end]),idxGroup,21)

    AligImg = abs(AligImg);
    save(['Rep8AligImg_',sMethod,'.mat'],'AligImg','-v7.3');
end