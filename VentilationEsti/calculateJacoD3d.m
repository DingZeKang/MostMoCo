function JacoBianDeter = calculateJacoD3d(motionField)

siz = size(motionField);
JacoBianDeter = zeros(siz(1:3));

% Grad11 = motionField([2:end,end-1],:,:,1)-motionField(:,:,:,1);
% Grad12 = motionField(:,[2:end,end-1],:,1)-motionField(:,:,:,1);
% Grad13 = motionField(:,:,[2:end,end-1],1)-motionField(:,:,:,1);
% 
% Grad21 = motionField([2:end,end-1],:,:,2)-motionField(:,:,:,2);
% Grad22 = motionField(:,[2:end,end-1],:,2)-motionField(:,:,:,2);
% Grad23 = motionField(:,:,[2:end,end-1],2)-motionField(:,:,:,2);
% 
% Grad31 = motionField([2:end,end-1],:,:,3)-motionField(:,:,:,3);
% Grad32 = motionField(:,[2:end,end-1],:,3)-motionField(:,:,:,3);
% Grad33 = motionField(:,:,[2:end,end-1],3)-motionField(:,:,:,3);


% [Grad11,Grad12,Grad13] = gradient(motionField(:,:,:,1));%X
% [Grad21,Grad22,Grad23] = gradient(motionField(:,:,:,2));%Y 
% [Grad31,Grad32,Grad33] = gradient(motionField(:,:,:,3));%Z
% 
% % Grad12 = Grad12+1;
% Grad11 = Grad11+1;
% % Grad13 = Grad13+1;
% Grad22 = Grad22+1;
% % Grad21 = Grad21+1;
% % Grad23 = Grad23+1;
% % Grad32 = Grad32+1;
% % Grad31 = Grad31+1;
% Grad33 = Grad33+1;
% 
% for i = 1:siz(1)
%    for j = 1:siz(2)
%       for k = 1:siz(3)
%          JacoMatrix = [Grad11(i,j,k),Grad12(i,j,k),Grad13(i,j,k);Grad21(i,j,k),Grad22(i,j,k),Grad23(i,j,k);...
%              Grad31(i,j,k),Grad32(i,j,k),Grad33(i,j,k)];
%          JacoBianDeter(i,j,k) =  det(JacoMatrix);
%       end
%    end
% end


% [Grad22,Grad21,Grad23] = gradient(motionField(:,:,:,2));%Y 1
% [Grad12,Grad11,Grad13] = gradient(motionField(:,:,:,1));%X 2
% [Grad32,Grad31,Grad33] = gradient(motionField(:,:,:,3));%Z 3

[Grad12,Grad11,Grad13] = gradient(motionField(:,:,:,2));%Y 1
[Grad22,Grad21,Grad23] = gradient(motionField(:,:,:,1));%X 2
[Grad32,Grad31,Grad33] = gradient(motionField(:,:,:,3));%Z 3, this is right manner, 20211122, DingZekang

% Grad12 = Grad12+1;
Grad11 = Grad11+1;
% Grad13 = Grad13+1;
Grad22 = Grad22+1;
% Grad21 = Grad21+1;
% Grad23 = Grad23+1;
% Grad32 = Grad32+1;
% Grad31 = Grad31+1;
Grad33 = Grad33+1;

for i = 1:siz(1)
   for j = 1:siz(2)
      for k = 1:siz(3)
         JacoMatrix = [Grad11(i,j,k),Grad12(i,j,k),Grad13(i,j,k);Grad21(i,j,k),Grad22(i,j,k),Grad23(i,j,k);...
             Grad31(i,j,k),Grad32(i,j,k),Grad33(i,j,k)];
         JacoBianDeter(i,j,k) =  det(JacoMatrix);
      end
   end
end