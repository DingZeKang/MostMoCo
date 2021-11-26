function TrueIndex = RespirationEsti_KC0(kc0,ThrowSeg,seqParam)
ntviews=size(kc0,1);
LCH =size(kc0,2);
spk = seqParam.spk;
seg = seqParam.seg;
TR = seqParam.TR;% unit ms
TrajDistr = seqParam.TrajDistr;
LPE = spk*seg;
ThrowLine=ThrowSeg*spk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Respiratory motion detection

t_N = size(kc0,1);

if TrajDistr ==3
    win_len = ceil(100/TR);
elseif TrajDistr ==4
    win_len = ceil(100/TR);
end
ksp_r = reshape(kc0(1:floor(t_N/win_len)*win_len,:),win_len,floor(t_N/win_len),size(kc0,2));
ksp_r = squeeze(sum(ksp_r,1));

[U,S,V ] = svd(ksp_r,'econ');
eig_num = 4;
eig_nums = 2;
a = U(:,eig_nums:eig_num)*diag(S(eig_nums:eig_num,eig_nums:eig_num));
a = mean(a(:));
a = a/abs(a);
if abs(real(a))>=abs(imag(a))
    Res_Signal = real(U(:,1:eig_num)*diag(S(1:eig_num,1:eig_num))*conj(a));
else
    Res_Signal = imag(U(:,1:eig_num)*diag(S(1:eig_num,1:eig_num))*conj(a));
end
Res_Signal = interp1(Res_Signal([1,1:end,end]),(1:t_N)'/win_len+1);
Res_Signal = filt1d(Res_Signal,TR/1000,[0.05,1],100);
figure,plot(Res_Signal(1:15000));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort the k-space data and trajectory according to respiratory motion 
[~,index]=sort(Res_Signal,'descend');
TrueIndex = index;
end

function k0f = filt1d(k0,TR,f0,N)
% 1d low pass filter
% Inputs:
%   k0: signal
%   TR: repetition time(s)
%   f0: motion freq(Hz)
%   N:  filter length
% Outputs:
%   k0f:filtered k0

if nargin<3
    f0 = 1;
end
if nargin<4
    N = 100;
end
N = ceil(max(2/(f0(1)*TR),N)/2)*2;
win = fir1(N,2*f0*TR,'bandpass');
k0e = zeros(length(k0(:))+N+1,1);
k0e(N/2+2:end-N/2) = k0(:);
k0e(1:N/2+1)=k0(1);
k0e(end-N/2+1:end)=k0(end);
k0f = conv(k0e(:),win,'same');
k0f = k0f(N/2+2:end-N/2);
end

function k0f = filt2d(k0,TR,f0,N)
% 1d high pass filter
% Inputs:
%   k0: signal
%   TR: repetition time(s)
%   f0: motion freq(Hz)
%   N:  filter length
% Outputs:
%   k0f:filtered k0

if nargin<3
    f0 = 1;
end
if nargin<4
    N = 100;
end
N = ceil(max(2/(f0*TR),N)/2)*2;
win = fir1(N,2*f0*TR,'high');
k0e = zeros(length(k0(:))+N+1,1);
k0e(N/2+2:end-N/2) = k0(:);
k0e(1:N/2+1)=k0(1);
k0e(end-N/2+1:end)=k0(end);
k0f = conv(k0e(:),win,'same');
k0f = k0f(N/2+2:end-N/2);
end