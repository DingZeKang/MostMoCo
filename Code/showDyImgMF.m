function showDyImgMF(img,idxGroup,figidx)
% Function to show motion field

traidx = idxGroup(1);
coridx = idxGroup(2);
sagidx = idxGroup(3);

% img = abs(img/max(abs(img(:))));
M1 = [];M2 = [];M3 = [];
img_sz = size(img);
if length(img_sz) ==3
   img_sz(4)=1; 
end
range1 = [];
range2 = [];
range3 = [];
for i = 1:img_sz(4)
    if traidx
        M1 = [M1,flipud(squeeze(img(:,:,traidx,i))')];
%         M1 = M1/max(abs(M1(:)));
        
    end
    
    
    if coridx
        M2 = [M2,fliplr((squeeze(img(:,coridx,:,i)))')];
%         M2 = M2/max(abs(M2(:)));
        
    end
    
    
    if sagidx
        M3 = [M3,fliplr((squeeze(img(sagidx,:,:,i)))')];
%         M3 = M3/max(abs(M3(:)));
        
    end
   
end
figure(figidx),subplot(311),imshow(M1,range1),colormap Jet,colorbar;
figure(figidx),subplot(312),imshow(M2,range2),colormap Jet,colorbar;
figure(figidx),subplot(313),imshow(M3,range3),colormap Jet,colorbar;