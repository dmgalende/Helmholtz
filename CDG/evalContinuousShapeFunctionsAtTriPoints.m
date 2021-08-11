function [shapeFunctions,derShapeFunctions] = ...
    evalContinuousShapeFunctionsAtTriPoints(x,nDeg,coord,op)

% given nPoints in the referenceTriangle evaluate the continuous shape
% functions of degree less or equal to nDeg at these points
%
% Inputs:
%   x :    points in (xi,eta) coordinates in the reference triangles
%   nDeg:  degree of interpolation
%   coord: nodal coordinates of the reference element
%
% Outputs:
%   shapeFunctions:    continuous shape functions evaluate at x
%                      size nDof X nOfPoints


N = (nDeg+1)*(nDeg+2)/2;
nOfPoints = size(x,1);
shapeFunctions = zeros(N,nOfPoints);
derShapeFunctions = zeros(N,nOfPoints,2);

V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

for i=1:nOfPoints
    [p,dp_dxi,dp_deta] = orthopoly2D_deriv_xieta(x(i,:),nDeg);
    shapeFunctions(:,i) = (U\(L\(P*p)))';
    derShapeFunctions(:,i,1) = (U\(L\(P*dp_dxi)))';
    derShapeFunctions(:,i,2) = (U\(L\(P*dp_deta)))';
end