function Fu = nurbsCurveDerivPoint(nurbs,u)
%
% Fu = nurbsCurveDerivPoint(nurbs,u)
%
% Only first derivatives at this moment
% (for higher derivatives see The NURBS Book - Page 125)
%
% Input:
% nurbs:  struct containing the nurbs curve information
% u:      parameter 
%
% Output:
% Fu:     derivative of the nurbs at point NURBS(u)
%

nurbsDer = nurbsCurveDerivControlPoints(nurbs);
[p1,w1] = nurbsCurvePoint(nurbsDer,u,'noHom');
[p2,w2] = nurbsCurvePoint(nurbs,u);
Fu = (p1 - p2*w1)/w2;