function uIni = initialGuessProjection(nurbs, xNode)

x = xNode(1);
y = xNode(2);

nParams = 10000;
us = linspace(nurbs.iniParam, nurbs.endParam, nParams);
Nus = zeros(nParams,2);
for i = 1:nParams
    pt = nurbsCurvePoint(nurbs, us(i));
    Nus(i,:) = pt(1:2);
end
dists = sqrt( (x - Nus(:,1)).^2 + (y - Nus(:,2)).^2  );
[m,p]=min(dists);
uIni = us(p);