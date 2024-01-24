function res = mtimes(a,b)

if a.adjoint
    % could design weights for each bin
    % bx = repmat(mean(b,5),1,1,1,1,size(b,5));
    bx = zeros(size(b));
    w_s = zeros(size(b,4),1);
    for i = 1:size(b,4)
        tmp = [1:size(b,4)]';
        w = (a.coef).^(-abs(tmp-i));
        w(i) = 0;
        w_s(i) = sum(w);
    end
    for i = 1:size(b,4)
        tmp = [1:size(b,4)]';
        w = (a.coef).^(-abs(tmp-i));
        w(i) = 0;
        w = w./w_s;
        w = repmat(permute(w,[2,3,4,1]),size(b,1),size(b,2),size(b,3),1);
        bt = movepixel(b,a.BAg(:,:,:,:,:,i));
%         showDyImg(squeeze(bt),[100 70 135],11)
        bt = sum((bt.*w),4);
        bx(:,:,:,i) = bt;
    end
    res = b - bx;
    
else
    % could design weights for each bin
    % bx = repmat(mean(b,5),1,1,1,1,size(b,5));
    bx = zeros(size(b));
    for i = 1:size(b,4)
        tmp = [1:size(b,4)]';
        w = (a.coef).^(-abs(tmp-i));
        w(i) = 0;
%         w = (abs(tmp-i)+1).^(-1);
        w_s = sum(w);
        w = repmat(permute(w,[2,3,4,1]),size(b,1),size(b,2),size(b,3),1);
        
        bt = movepixel(b,a.BAg(:,:,:,:,:,i));%per motion state -> ith motion state
%         showDyImg(squeeze(bt),[100 70 135],11)
        bt = sum((bt.*w/w_s),4);
        bx(:,:,:,i) = bt;
    end
    res = b - bx;
%     showDyImg(res(:,:,:,[1,5]),[250 95 160],900)
end