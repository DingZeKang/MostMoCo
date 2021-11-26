function y = Forward(x,param,E,reg_field,weighting)
% DingZekang SJTU
% 2020 07 at UIH
y = zeros(size(param.y));
y = permute(y,[5,1,2,3,4]);
dyImage = zeros([size(x),size(y,5)]);

for i = 1:size(y,5)
    dyImage(:,:,:,i) = imwarp(x,-reg_field(:,:,:,:,i));
    ksp = E{i}*dyImage(:,:,:,i);
    ksp = ksp * weighting(i);
    y(1,:,:,:,i) = reshape(ksp,[size(param.y,1),size(param.y,2),size(param.y,3)]);
end
