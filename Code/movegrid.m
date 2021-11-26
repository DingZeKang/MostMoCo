function Fxy = movegrid(F1x,B)
% Forward x-y based on x grid, x deformed grid, y target grid
% F: F1->x
% B: Bx->1 - By->1
Fxy = zeros(size(B));
for i = 1:size(B,5)
    for j = 1:size(B,4)
        Fxy(:,:,:,j,i) =  imwarp(B(:,:,:,j,i),F1x(:,:,:,:,i));
    end
end