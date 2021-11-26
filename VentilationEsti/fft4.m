function ksp = fft4(image)
siz = size(image);
ksp = zeros(siz);
for i = 1:siz(4)
   ksp(:,:,:,i) = fftshift(fftn(squeeze(image(:,:,:,i)))); 
end