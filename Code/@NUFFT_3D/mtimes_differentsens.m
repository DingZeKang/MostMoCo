function ress = mtimes(a,bb)

 if a.adjoint        % iNUFFT
     % Multicoil non-Cartesian k-space to Cartesian image domain
     % nufft for each coil and time point
     res = zeros([size(a.b1),size(bb,4)]);
     for tt=1:size(bb,4)
        for ch=1:size(a.b1,4)
            b = bb(:,:,ch,tt) .* a.w(:,:,tt);
            %res(:,:,:,tt) = reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize)),a.imSize(1),a.imSize(2),a.imSize(3));%origin
            res(:,:,:,ch,tt) = reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize)),a.imSize(1),a.imSize(2),a.imSize(3));
        end
     end
     % compensate for undersampling factor
%      res=res*size(a.b1,1)*pi/2/size(a.w,2);     
     % coil combination for each time point
     ress = zeros(size(a.b1,1),size(a.b1,2),size(a.b1,3),size(bb,4));
     for tt=1:size(bb,4)
%          ress(:,:,:,tt)=sum(res(:,:,:,tt).*conj(a.b1),4)./sum(abs((a.b1)).^2,4); %#ok<AGROW>%%origin
         ress(:,:,:,tt)=sum(res(:,:,:,:,tt).*conj(a.b1(:,:,:,:,tt)),4)./sum(abs((a.b1(:,:,:,:,tt))).^2,4);
     end
 else
     % Cartesian image to multicoil non-Cartesian k-space   NUFFT 
     ress = zeros([a.dataSize,size(a.b1,4),size(bb,4)]);
     for tt=1:size(bb,4)
        for ch=1:size(a.b1,4)
            res=bb(:,:,:,tt).*a.b1(:,:,:,ch,tt); 
            %ress(:,:,ch,tt) = reshape(nufft(res,a.st{tt})/sqrt(prod(a.imSize)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,tt);%origin
            ress(:,:,ch,tt) = reshape(nufft(res,a.st{tt})/sqrt(prod(a.imSize)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,tt);
        end
     end
 end