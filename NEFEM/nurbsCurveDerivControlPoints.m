function nurbsDer = nurbsCurveDerivControlPoints(nurbs) 
%
% nurbsDer = nurbsCurveDerivControlPoints(nurbs) 
%

U = nurbs.U;
Pw = nurbs.Pw;
p = length(find(U==U(1))) - 1;
m = length(U) - 1;
n = m - p - 1;

nurbsDer.U = U(2:end-1);

nurbsDer.Pw = zeros(n,4);
for i = 1:n
    nurbsDer.Pw(i,:) = p*(Pw(i+1,:)-Pw(i,:))/(U(i+p+1)-U(i+1));    
end