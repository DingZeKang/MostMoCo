function kc0 = Extractk0(FIDrawdata,seqParam)

[RO,LPE,nc]=size(FIDrawdata);
spk = seqParam.spk;
UTEDataBuffer = seqParam.UTEDataBuffer;
TrajDistr  = seqParam.TrajDistr;
if TrajDistr == 4
%     ntviews = LPE/spk;
    ntviews = LPE;
else% 0 3
    ntviews = LPE;
end
kc0 = zeros(ntviews,nc);

kc0 = squeeze(FIDrawdata(1:UTEDataBuffer+1,:,:));
kc0 = squeeze(mean(kc0,1));