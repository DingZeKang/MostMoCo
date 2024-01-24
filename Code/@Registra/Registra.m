function res =Registra(L,B,F)
%  3DTV motion constraint

res.F = F;
res.B = B;
res.adjoint = 0;
res.w = ones(L,1);
res = class(res,'Registra');