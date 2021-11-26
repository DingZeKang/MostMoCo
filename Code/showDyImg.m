function showDyImg(img,idxGroup,figidx)

traidx = idxGroup(1);
coridx = idxGroup(2);
sagidx = idxGroup(3);

img = double(img);
img = abs(img/max(abs(img(:))));
M1 = [];M2 = [];M3 = [];
img_sz = size(img);
if length(img_sz) ==3
   img_sz(4)=1; 
end
range = [0,0.25];
for i = 1:img_sz(4)
    if traidx
        M1 = [M1,flipud(squeeze(img(:,:,traidx,i))')];
%         M1 = M1/max(abs(M1(:)));
        figure(figidx),subplot(311),imshow(M1,range);
    end
    
    if coridx
        M2 = [M2,fliplr((squeeze(img(:,coridx,:,i)))')];
%         M2 = M2/max(abs(M2(:)));
        figure(figidx),subplot(312),imshow(M2,range);
    end
    
    if sagidx
        M3 = [M3,fliplr((squeeze(img(sagidx,:,:,i)))')];
%         M3 = M3/max(abs(M3(:)));
        figure(figidx),subplot(313),imshow(M3,range);
    end
end