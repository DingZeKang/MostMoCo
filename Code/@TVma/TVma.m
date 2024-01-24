function res =TVma(L,B,F,coef)
%  3DTV motion constraint

res.F = F;%ref->i
res.B = B;%i->ref
res.adjoint = 0;
res.nbin = L;%L is nbin
res.coef = coef;

res.BAg = B;
% res.BAg = zeros([size(B),L]);
% for i = 1:L
%     Bz = B-repmat(B(:,:,:,:,i),1,1,1,1,size(B,5));
%     F = repmat(F(:,:,:,:,i),1,1,1,1,size(F,5));
%     Fxy = movegrid(F,Bz);
%     res.BAg(:,:,:,:,:,i) = Fxy;
% end


res = class(res,'TVma');

