function res = mtimes(a,b)

if a.adjoint
%     res = adjDx(b(:,:,:,:,1)) + adjDy(b(:,:,:,:,2)) + adjDz(b(:,:,:,:,3)) + adjDt(b(:,:,:,:,4));
    res = adjDt(b);
else
%     Dx = b([2:end,end],:,:,:) - b;
%     Dy = b(:,[2:end,end],:,:) - b;
%     Dz = b(:,:,[2:end,end],:) - b;
%     Dt = b(:,:,:,[2:end,end]) - b;
%     res = cat(5,Dx,Dy,Dz,Dt); % cat: connect matrix at the 5th dimension
    res=b(:,:,:,[2:end,end])-b;
end

% function res = adjDx(x)
% res = x([1,1:end-1],:,:,:) - x;
% res(1,:,:,:) = -x(1,:,:,:);
% res(end,:,:,:) = x(end-1,:,:,:);
% 
% function res = adjDy(x)
% res = x(:,[1,1:end-1],:,:) - x;
% res(:,1,:,:) = -x(:,1,:,:);
% res(:,end,:,:) = x(:,end-1,:,:);
% 
% function res = adjDz(x)
% res = x(:,:,[1,1:end-1],:) - x;
% res(:,:,1,:) = -x(:,:,1,:);
% res(:,:,end,:) = x(:,:,end-1,:);

function res = adjDt(x)
res = x(:,:,:,[1,1:end-1]) - x;
res(:,:,:,1) = -x(:,:,:,1);
res(:,:,:,end) = x(:,:,:,end-1);