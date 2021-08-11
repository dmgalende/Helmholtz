function [pt,w] = nurbsCurvePoint(nurbs,u,strHom)
%
% [pt,w] = nurbsCurvePoint(nurbs,u,strHom)
% 
% Input:
% nurbs:  struct containing the nurbs curve information
% u:      parameter 
% strHom: (optional) if it's specified homogeneous 
%         coordinates are considered
%
% Output:
% pt:     point of the nurbs curve
% w:      weight (non-homogeneous coordinates)

U = nurbs.U;
Pw = nurbs.Pw;
p = length(find(U==U(1))) - 1;

span = nurbsCurveFindSpan(u,U);
N = nurbsCurveBasisFuns(span,u,U);

pt = 0;
for i=1:p+1
    pt = pt + N(i)*Pw(span-p+i-1,:);
end

w = pt(4);
if nargin==2
    pt = pt(1:3)/w;
else
    pt = pt(1:3);
end