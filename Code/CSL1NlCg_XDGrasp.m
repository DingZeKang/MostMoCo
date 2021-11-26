function x = CSL1NlCg_XDGrasp(x0,param)
% f(x) = ||E*x - y||^2 + L1Weight * ||W*x||_1TVWeight * ||TV*x||_1 + TV_TempWeight * ||TV_Temp*x||_1 

% starting point
x = x0;

% line search parameters
maxlsiter = 10 ;
gradToll = 1e-8 ;
param.l1Smooth = 1e-15;	
alpha = 0.01;  
beta = 0.6;
t0 = 1 ; 
k = 0;
% compute g0  = grad(f(x))
g0 = grad(x,param);
dx = -g0;

% iterations
while(1)

    % backtracking line-search
	f0 = objective(x,dx,0,param);
	t = t0;
    f1 = objective(x,dx,t,param);
	lsiter = 0;
    
    %line search
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:))) && (lsiter < maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(x,dx,t,param);
	end

% 	if lsiter == maxlsiter
% 		disp('Error - line search ...');
% 		return;
% 	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta; end 
	if lsiter < 1, t0 = t0 / beta; end

    % update x
	x = (x + t*dx);
    
	% print some numbers for debug purposes	
    if param.display
        fprintf('%d   , obj: %f, L-S: %d\n', k,f1,lsiter);
    end
    k = k + 1;
	
    %imshow
    if 0 
        p=round(110*2/3);
        figure, imshow(fliplr(abs(squeeze(x(:,p,:,8)))'),[])
    end
    
	% stopping criteria (to be improved)
	if (k >= param.nite) || (norm(dx(:)) < gradToll), break; end

    
    %conjugate gradient calculation
	g1 = grad(x,param);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk*dx;
	

end
return;

function res = objective(x,dx,t,param) %**********************************
m = x+t*dx;
% L2-norm part
w = param.F (m) - param.y;
L2Obj = w(:)' * w(:);

% TV part along time
if param.TV_TempWeightT
    w = param.TV_TempT * (m); 
    TV_TempObjT = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TV_TempObjT = 0;
end    

% TV part along spatial
if param.TV_TempWeightS
    w = param.TV_TempS * (m); 
    TV_TempObjS = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TV_TempObjS = 0;
end  

res = L2Obj + param.TV_TempWeightT*TV_TempObjT + param.TV_TempWeightS*TV_TempObjS;


function g = grad(x,param)%***********************************************
m=x;
% L2-norm part
L2Grad = 2.*( param.Ft ( param.F (m)-param.y ) );

% TV part along time
if param.TV_TempWeightT
    w = param.TV_TempT * m;
    TV_TempGradT = param.TV_TempT'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TV_TempGradT = 0;
end

% TV part along spatial
if param.TV_TempWeightS
    w = param.TV_TempS * m;
    TV_TempGradS = param.TV_TempS'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TV_TempGradS = 0;
end

g = L2Grad + param.TV_TempWeightT*TV_TempGradT + param.TV_TempWeightS*TV_TempGradS;