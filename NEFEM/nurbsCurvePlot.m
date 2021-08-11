function h = nurbsCurvePlot(nurbs, iniParam, endParam, nPoints)
%
% nurbsCurvePlot(nurbs, iniParam, endParam, nPoints)
%
% Input:
% nurbs: struct containing the nurbs curve information
% iniParam, endParam:   trimmed NURBS  (complete NURBS by default)
% nPoints:              Number of points (200 by default)
%

% By default options
if nargin < 2
    iniParam = nurbs.iniParam;
    endParam = nurbs.endParam;
end

if nargin < 4
    nPoints = 5000;
end

% Curve points
fu = zeros(nPoints,3);
uVector = linspace(iniParam,endParam,nPoints);
for iPoint = 1:nPoints
    fu(iPoint,:) = nurbsCurvePoint(nurbs,uVector(iPoint));
end

h = plot(fu(:,1), fu(:,2), 'k-', 'LineWidth', 1);