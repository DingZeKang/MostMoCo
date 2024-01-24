function res =TVs(dim)
% improved TV constraint

res.adjoint = 0;
res.dim = dim;
res = class(res,'TVs');