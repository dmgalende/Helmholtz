function  [xi,J] = inverseLinearMapping(Xv,x)
%
% [xi,J] = inverseLinearMapping(Xv,x)
% 

nOfSpatialDimensions = size(Xv,2);
switch nOfSpatialDimensions
    case 2
        [xi,J] = inverseLinearMapping_2D(Xv,x);
    case 3
        [xi,J] = inverseLinearMapping_3D(Xv,x);
    otherwise
        error('inverseLinearMapping: wrong nOfSpatialDimensions')
end
    
function [xi,J] = inverseLinearMapping_2D(Xv,x)

J = [Xv(2,1)-Xv(1,1) Xv(3,1)-Xv(1,1); 
     Xv(2,2)-Xv(1,2) Xv(3,2)-Xv(1,2)];

% Point in local coordinates
xi = J\[x(:,1)'-Xv(1,1); 
        x(:,2)'-Xv(1,2)];
xi = 2*xi'-1;

function [xi,J] = inverseLinearMapping_3D(Xv,x)

J = [Xv(2,1)-Xv(1,1), Xv(3,1)-Xv(1,1), Xv(4,1)-Xv(1,1) 
     Xv(2,2)-Xv(1,2), Xv(3,2)-Xv(1,2), Xv(4,2)-Xv(1,2) 
     Xv(2,3)-Xv(1,3), Xv(3,3)-Xv(1,3), Xv(4,3)-Xv(1,3)];
 
b = [x(:,1)'-Xv(1,1);
     x(:,2)'-Xv(1,2);
     x(:,3)'-Xv(1,3)];

 % Point in local coordinates
xi = J\b;
xi = 2*xi'-1;