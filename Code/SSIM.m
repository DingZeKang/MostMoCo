function [ssim,stdSSIM] = SSIM(ref, X)
% Function SSIM is to calculate SSIM (Structure Similarity)for qualitatively evaluating the image quality.
%
% Input Variables:
% ref - reference image , or called the ground truth 
% X - image to be evaluated
%
% Output Variables:
% ssim - Structure Similarity
% stdSSIM - standard deviation of SSIM (only for 3D case)
%
% Record of Revisions:
% July-31-2020===Zhao He===original code
% 04-23 2021 DingZekang Add 4D


n = ndims(ref); 
ref = abs(ref)/max(abs(ref(:)))*255;
X = abs(X)/max(abs(X(:)))*255;

if n==2 % if 2d data
    ssim = SSIM_2D(ref, X);
    stdSSIM = 0;
elseif n==3 % if 3d data
    for ii = 1:size(ref,3), E(ii)= SSIM_2D(ref(:,:,ii), X(:,:,ii)); end
    ssim = mean(E);  stdSSIM = std(E); 
elseif n == 4
   for ii = 1:size(ref,3)
      for jj = 1:size(ref,4)
          E(ii,jj)= SSIM_2D(ref(:,:,ii,jj), X(:,:,ii,jj));
      end
   end
   ssim = mean(E(:));  stdSSIM = std(E(:)); 
else
     disp('Input must be 2D ,3D or 4D data!');return;
end

end

function [mssim, ssim_map,siga_sq,sigb_sq] = SSIM_2D(ima, imb)
% ========================================================================
%ssim���㷨��Ҫ�ο��������ģ�
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
% quality assessment: From error visibility to structural similarity,"
% IEEE Transactios on Image Processing, vol. 13, no. 4, pp. 600-612,
% Apr. 2004.
%  ���ȶ�ͼ��Ӵ�����w=fspecial('gaussian', 11, 1.5);
%                 (2*ua*ub+C1)*(2*sigmaa*sigmab+C2)
%   SSIM(A,B)=������������������������������������������������
%              (ua*ua+ub*ub+C1)(sigmaa*sigmaa+sigmab*sigmab+C2)
%     C1=��K1*L��;
%     C2=(K2*L);   K1=0.01��K2=0.03
%     LΪ�Ҷȼ�����L=255
%-------------------------------------------------------------------
%     ima - �Ƚ�ͼ��A
%     imb - �Ƚ�ͼ��B
%
% ssim_map - ���Ӵ���õ���SSIM��A,B|w����ɵ�ӳ�����
%    mssim - �ԼӴ��õ���SSIM��A,B|w����ƽ���������յ�SSIM��A,B��
%  siga_sq - ͼ��A�������ڻҶ�ֵ�ķ���
%  sigb_sq - ͼ��B�������ڻҶ�ֵ�ķ���
%-------------------------------------------------------------------
%  Cool_ben
%========================================================================

w = fspecial('gaussian', 11, 1.5);	%window �Ӵ�
K(1) = 0.01;					
K(2) = 0.03;					
L = 255;     
ima = double(ima);
imb = double(imb);

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
w = w/sum(sum(w));

ua   = filter2(w, ima, 'valid');%�Դ����ڲ�û�н���ƽ�������������˹�����
ub   = filter2(w, imb, 'valid'); % ���Ƽ�Ȩƽ��
ua_sq = ua.*ua;
ub_sq = ub.*ub;
ua_ub = ua.*ub;
siga_sq = filter2(w, ima.*ima, 'valid') - ua_sq;
sigb_sq = filter2(w, imb.*imb, 'valid') - ub_sq;
sigab = filter2(w, ima.*imb, 'valid') - ua_ub;

ssim_map = ((2*ua_ub + C1).*(2*sigab + C2))./((ua_sq + ub_sq + C1).*(siga_sq + sigb_sq + C2));


mssim = mean2(ssim_map);

return
end