function [B,F] = registrationPandG(ref,nbin,Img,imgsize,UsingGroupReg,seqParam)

Img = abs(imgauss4d(Img,.5));
[Img] = normalDifSta(Img);

B = zeros([size(Img(:,:,:,1)),3,nbin,nbin]);
F = zeros([size(Img(:,:,:,1)),3,nbin]);
if ~ UsingGroupReg
    traidx = 250;
    coridx = 95;
    sagidx = 160;
    idxGroup = round([traidx,coridx,sagidx]);
    Img = double(Img);
    Img = Img/max(Img(:));
    Imgeq = zeros(size(Img));
    ImgReg = zeros(size(Img));
    for ref = 1:nbin
        fixedGPU = gpuArray(Img(:,:,:,ref));
        fixedHist = imhist(fixedGPU);
        for i = 1:nbin

            movingGPU = gpuArray(Img(:,:,:,i));   
            movingGPU = histeq(movingGPU,fixedHist);
            [B_GPU,movingRegGPU] = imregdemons(movingGPU,fixedGPU,[200,100,50,25]*4,'PyramidLevels',4,'DisplayWaitbar',false,'AccumulatedFieldSmoothing',1.0);
            B(:,:,:,:,i,ref) = gather(B_GPU);

        end
    end

else
    
    % ToDo, GroupWise Registration
    
end