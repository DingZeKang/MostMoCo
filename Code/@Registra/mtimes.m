function res = mtimes(a,b)

if a.adjoint
    res = zeros(size(b));
    % could design weights for each bin
    % bx = repmat(mean(b,5),1,1,1,1,size(b,5));
    M_state = 4;
    for i = 1:size(b,M_state)
        res(:,:,:,i) = imwarp(b(:,:,:,i),a.F(:,:,:,:,i));
    end
else
    res = zeros(size(b));
    % could design weights for each bin
    % bx = repmat(mean(b,5),1,1,1,1,size(b,5));
    M_state = 4;
    for i = 1:size(b,M_state)
        res(:,:,:,i) = imwarp(b(:,:,:,i),a.B(:,:,:,:,i));
    end
end