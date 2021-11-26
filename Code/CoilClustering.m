function [dave, coilID] = CoilClustering(d1, thresh)
 
% Coil clustering
%
% Input:
% (1) d1, navigator signal 
% (2) thresh, correlation threshold
%
% Output:
% (1) dave, averaged navigator within the cluster
% (2) coilID, mask indicating whether the coil element is selected in the cluster
%
% (c) Tao Zhang, 2015, Stanford University

[nviews, nc] = size(d1);

% find out the covariance matrix between coils
% disp('Calculating covariance matrix...');
corrm = zeros(nc,nc);

for i = 1 : nc
	for j = i : nc
		corrm(i,j) = xcov(d1(:,i),d1(:,j),0,'coef');
		corrm(j,i) = corrm(i,j);
	end
end

% disp('Spectral clustering...');
% set a mask according to this value
mask = zeros(nc,nc);
mask(abs(corrm)>thresh) = 1;

[U,S,V] = svd(mask,'econ');
v1 = abs(U(:,1));

thresh2 = 0.1;
subgroup = zeros(nc,1);
subgroup(v1>thresh2) = 1;

dave = zeros(nviews,1);

% adjust from the first coil
subindex = find(subgroup==1);
coilID = subgroup;

if( sum(subgroup) < 2)
    disp('Clustering failed, only one coil was selected...');
end

for c = 1:nc
	if(subgroup(c)>0)
		if(corrm(subindex(1),c)>0)
			dave = dave + d1(:,c);
		else
			dave = dave - d1(:,c);
			coilID(c) = -coilID(c);
		end
	end
end

dave = dave./sum(subgroup);


