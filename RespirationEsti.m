function [TrueIndex,Res_Signal] = RespirationEsti(kdata_SI,ThrowSeg,seqParam)
ntviews=size(kdata_SI,2) - ThrowSeg;
LCH =size(kdata_SI,3);
spk = seqParam.spk;
seg = seqParam.seg;
LPE = spk*seg;
ThrowLine=ThrowSeg*spk;
TR = seqParam.TR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Respiratory motion detection

kc = kdata_SI(1:end,:,:);

ZIP=abs((ifft((kc),size(kc,1),1)));  

for ii=1:LCH
    for jj=1:ntviews
        maxprof=max(ZIP(:,jj,ii));
        minprof=min(ZIP(:,jj,ii));
        ZIP(:,jj,ii)=(ZIP(:,jj,ii)-minprof)./(maxprof-minprof);
    end
end
figure,imagesc(abs(ZIP(:,:,5))),axis image,colormap(gray), axis off

kk=1;clear PCs
for ii=1:LCH
    tmp=permute(ZIP(:,:,ii),[1,3,2]);
    tmp=abs(reshape(tmp,[size(tmp,1)*size(tmp,2),ntviews])');
    
    covariance=cov(tmp);
    covariance(isnan(covariance)) = 0;
    [tmp2, V] = eig(covariance);
    V = diag(V);
    [junk, rindices] = sort(-1*V);
    V = V(rindices);
    tmp2 = tmp2(:,rindices);
    PC = (tmp2' * tmp')';
    % Take the first two principal components from each coil element.
    for jj=1:1
        
        tmp3=smooth(PC(:,jj),6,'lowess'); % do some moving average smoothing
        
        tmp3=tmp3-min(tmp3(:));
        tmp3=tmp3./max(tmp3(:));       
        
        PCs(:,kk)=tmp3;kk=kk+1;
    end
end              
% close all
% Do coil clusting to find the respiratory motion signal
% Function obtained from Tao Zhang (http://web.stanford.edu/~tzhang08/software.html)
% figure,plot(-PCs(1:400,16));
thresh = 0.95;%origin 0.95
[Res_Signal, cluster] = CoilClustering(PCs, thresh);

%Normalize the signal for display
Res_Signal=Res_Signal-min(Res_Signal(:));
Res_Signal=Res_Signal./max(Res_Signal(:));
figure,plot(Res_Signal(1:end))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort the k-space data and trajectory according to respiratory motion 
[~,index]=sort(Res_Signal,'descend');
%figure,plot(Res_Signal,'r');
TrueIndex=1:LPE-ThrowLine;
A=1:spk;
for i=1:seg-ThrowSeg
    TrueIndex((i-1)*spk+1:i*spk) =(index(i)-1)*spk+A;
end
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
win = fir1(N,2*f0*TR,'bandpass');%'high' 'low' 'bandpass' 'stop'
k0e = zeros(length(k0(:))+N+1,1);
k0e(N/2+2:end-N/2) = k0(:);
k0e(1:N/2+1)=k0(1);
k0e(end-N/2+1:end)=k0(end);
k0f = conv(k0e(:),win,'same');
k0f = k0f(N/2+2:end-N/2);
end