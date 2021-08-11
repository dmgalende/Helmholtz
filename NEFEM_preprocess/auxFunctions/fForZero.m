function H = fForZero(u, xNode, nurbs)

if u<nurbs.iniParam
    u = nurbs.iniParam;
elseif u>nurbs.endParam
    u = nurbs.endParam;
end

point = nurbsCurvePoint(nurbs, u);
pointD = nurbsCurveDerivPoint(nurbs, u);
point = point(1:2);
pointD = pointD(1:2);
H = (xNode-point)*pointD';
