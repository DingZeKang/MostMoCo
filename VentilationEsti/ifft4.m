function image = ifft4(ksp)
siz = size(ksp);
image = zeros(siz);
for i = 1:siz(4)
   image(:,:,:,i) = ifftn(ifftshift(squeeze(ksp(:,:,:,i)))); 
end