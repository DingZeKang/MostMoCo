function [DyLung1] = normalDifSta(DyLung)

DyLung1 = zeros(size(DyLung));

for i = 1:size(DyLung,4)
%     tmp = DyLung(:,:,:,i);
    tmp = DyLung(:,:,:,:);
    DyLung1(:,:,:,i) = abs(DyLung(:,:,:,i))/max(abs(tmp(:)));

end