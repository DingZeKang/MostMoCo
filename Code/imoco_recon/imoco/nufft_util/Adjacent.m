function x = Adjacent(y,E,reg_field,weighting)
% DingZekang SJTU
% 2020 07 at UIH
img_size = size(reg_field);
img_size = img_size([1,2,3,5]);
x = zeros(img_size);

[RO1,LPE,LCH_PCA,tt] = size(y);

y = reshape(y,[RO1*LPE,LCH_PCA,tt]);
for i = 1:tt
    tmp = E{i}'*(y(:,:,i));
    x(:,:,:,i) = imwarp(tmp,reg_field(:,:,:,:,i));
end
x = mean(x,4)*tt/sum(weighting);