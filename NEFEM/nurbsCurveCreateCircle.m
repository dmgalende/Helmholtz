function nurbs = nurbsCurveCreateCircle(center,r,iniParam,endParam)
%
% nurbs = nurbsCurveCreateCircle(center,r,iniParam,endParam)
%
% Input:
% center, r: center and radius
% iniParam, endParam: Parameters for the trimmed nurbs
% Default option [iniParam, endParam]=[0,1] 
% 
% Output:
% nurbs: struct containing the nurbs information
%

U = [0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1];
p = 2;
m = size(U,2) - 1;
n = m-p-1;
P =[-1  0 0
    -1 -1 0
     0 -1 0
     1 -1 0
     1  0 0
     1  1 0
     0  1 0
    -1  1 0
    -1  0 0];

P = P*r;
P(:,1) = P(:,1) + center(1);
P(:,2) = P(:,2) + center(2);

b = sqrt(2)/2;
w = [1,b,1,b,1,b,1,b,1]';

Pw = zeros(n+1,4);
Pw(:,1) = P(:,1).*w;
Pw(:,2) = P(:,2).*w;
Pw(:,3) = P(:,3).*w;
Pw(:,4) = w;
    

% Se guarda la info de la nurbs en un struct
 nurbs.U = U;
 nurbs.Pw = Pw;
 
 if nargin == 2 
     iniParam = 0;
     endParam = 1;
 end
 
 nurbs.iniParam = iniParam;
 nurbs.endParam = endParam;
     