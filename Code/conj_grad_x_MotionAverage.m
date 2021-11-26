function x_k1 = conj_grad_x_MotionAverage(A,At,x_k0,d,rho,mTVx,z_k1,u_k0,rho1,mTVs,z1_k1,u1_k0,tol,cg_iterM)

% DingZekang SJTU
% 2021 04, at Med-X, SJTU

vec= @(x) x(:).' ;%row vector
col= @(x) x(:);%column vector

ResX = @(x) reshape(x,size(x_k0));
Resd= @(x) reshape(x,size(d));
ResZ= @(x) reshape(x,size(z_k1));
ResZ1= @(x) reshape(x,size(z1_k1));

    
    
b=[col(d(:));col(sqrt(rho/2)*(z_k1-u_k0/rho));col(sqrt(rho1/2)*(z1_k1-u1_k0/rho1))];
init=col(x_k0);

x=lsqr(@afun2,b,tol,cg_iterM,[],[],init); %add initial guess
    function y= afun2(x,transp_flag)
        if strcmp(transp_flag,'notransp')
            y1=col(A(ResX(x)));
            y2=sqrt(rho/2).*col(mTVx*(ResX(x)));
            y3=sqrt(rho1/2).*col(mTVs*(ResX(x)));
            y=[y1;y2;y3];
        else
            x1=x(1:numel(d));
            x2=x((numel(d)+1):(numel(d)+numel(z_k1)));
            x3=x((numel(d)+numel(z_k1)+1):end);
            y1=(At(Resd(x1)));
            y2=sqrt(rho/2).*(mTVx'*(ResZ(x2)));
            y3=sqrt(rho1/2).*(mTVs'*(ResZ1(x3)));
            y=col(y1+y2+y3);
        end
    end

x_k1=ResX(x);

end