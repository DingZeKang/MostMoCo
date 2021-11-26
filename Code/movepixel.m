function Iout = movepixel(I,B)
Iout = zeros(size(I));
for i = 1:size(I,ndims(I))
    Iout(:,:,:,i) =  imwarp(I(:,:,:,i),B(:,:,:,:,i));
end