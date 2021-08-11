function ZI = planarInterpolation(X,Y,Z,XI,YI)
% planar interpolation of the function Z
% X: 3 x 1
% Y: 3 x 1
% Z: 3 x 1
aux = [3 1];
if any(size(X)~=aux) || any(size(Y)~=aux) || size(Z,1)~=3 ...
        || any(size(XI)~=size(YI))
    error('wrong input dimension')
end
A = [X,Y,ones(3,1)];
plane = A\Z;
ZI = zeros(size(XI,1),size(Z,2));
for i = 1:size(plane,2)
    ZI(:,i) = plane(1,i)*XI + plane(2,i)*YI + plane(3,i);
end