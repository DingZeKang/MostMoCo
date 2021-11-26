function [tu,wu,kdatau] = data_sorting(traj,w,ksdata,nu,resp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data sorting
[nx,ns,nc] = size(ksdata);

% nu:  number of respiratory phases
nline = floor(ns/nu);

% Sort the k-space data and trajectory according to respiratory motion
if nargin > 4
    [~,index] = sort(resp,'descend');
    ksdata = ksdata(:,index,:);
    traj = traj(:,:,index);  
    w = w(:,index);
end

tu = zeros(3,nx,nline,nu);
wu = zeros(nx,nline,nu);
kdatau = zeros(nx,nline,nc,nu);
for ii = 1 : nu
    tu(:,:,:,ii) = traj(:,:,(ii-1)*nline+1:ii*nline);  
    wu(:,:,ii) = w(:,(ii-1)*nline+1:ii*nline);
    kdatau(:,:,:,ii) = ksdata(:,(ii-1)*nline+1:ii*nline,:);
end
