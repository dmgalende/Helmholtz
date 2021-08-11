function J = Jacobian(X)
%
% J = Jacobian(X)
%
% Jacobian of a given simplex
%

nsd = size(X,2);
switch nsd
    case 1
        J = (X(2)-X(1))/2;
    case 2
        J = Jacobian2D(X);
    case 3
        J = Jacobian3D(X);    
end

function J = Jacobian2D(X)
v1 = X(1,:);  v2 = X(2,:);  v3 = X(3,:);
J = [(v2-v1)/2 ; (v3-v1)/2];

function J = Jacobian3D(X)
v1 = X(1,:);  v2 = X(2,:);  v3 = X(3,:);  v4 = X(4,:);
J = [(v2-v1)/2 ; (v3-v1)/2 ; (v4-v1)/2];


