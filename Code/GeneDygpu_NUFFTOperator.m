function Dygpu_NUFFTOperator = GeneDygpu_NUFFTOperator(k_cs,w_cs,sens1,sw)
%DingZekang 20200711 at UIH
% k_cs      3 RO LPE nt
% w_cs      RO LPE nt
[~,RO1,LPE,tt] = size(k_cs);

img_size = size(sens1);
img_size = img_size(1:3);

k_cs = reshape(k_cs,[3,RO1*LPE,tt]);
w_cs = reshape(w_cs,[RO1,LPE,tt]);

Dygpu_NUFFTOperator = cell(1,tt);
for i = 1:tt
    k = k_cs(:,:,i);
    w = w_cs(:,:,i);
%     max(abs(w(:)))
    Dygpu_NUFFTOperator{i} = gpuNUFFT(k,w,2,3,sw,img_size,sens1,true);
end