function res = mtimes(a,b)

%eidt by DingZekang SJTU

if a.adjoint
    res = adjDx(b(:,:,:,:,1)) + adjDy(b(:,:,:,:,2)) + adjDz(b(:,:,:,:,3));
else
    Dx = b([2:end,end],:,:,:) - b;
    Dy = b(:,[2:end,end],:,:) - b;
    Dz = b(:,:,[2:end,end],:) - b;
    res = cat(5,Dx,Dy,Dz); % cat: connect matrix at the 5th dimension
end

function res = adjDx(x)
res = x([1,1:end-1],:,:,:) - x;
res(1,:,:,:) = -x(1,:,:,:);
res(end,:,:,:) = x(end-1,:,:,:);

function res = adjDy(x)
res = x(:,[1,1:end-1],:,:) - x;
res(:,1,:,:) = -x(:,1,:,:);
res(:,end,:,:) = x(:,end-1,:,:);

function res = adjDz(x)
res = x(:,:,[1,1:end-1],:) - x;
res(:,:,1,:) = -x(:,:,1,:);
res(:,:,end,:) = x(:,:,end-1,:);