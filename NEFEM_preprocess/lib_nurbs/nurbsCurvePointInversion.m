function parameter = nurbsCurvePointInversion(xNode, nurbs)
%
% parameter = nurbsCurvePointInversion(xNode, nurbs)
%

x = xNode(1);
y = xNode(2);

% In order to choose a good initial approximation the parametric space is
% divided into nParams intervals. The minimum distance is used for
% selecting the initial guess.
nParams = 10000;

% For Newton-Raphson
global tol tolParam
nIterMax = 1000;

% Initial approximation
us = linspace(nurbs.iniParam, nurbs.endParam, nParams);
nParams = length(us);
Nus = zeros(nParams,2);
for i = 1:nParams
    pt = nurbsCurvePoint(nurbs, us(i));
    Nus(i,:) = pt(1:2);
end
dists = sqrt( (x - Nus(:,1)).^2 + (y - Nus(:,2)).^2  );
[m,p]=min(dists);

% Look a root of the distance function
u = us(p);
if m<tol
    parameter=u;
else
    j=1;
    while j<nIterMax
        % Corrections
        if u < nurbs.iniParam
            u = nurbs.iniParam;
        elseif u > nurbs.endParam
            u = nurbs.endParam;
        end
        
        % Point and derivative (for Newton-Raphson)
        N = nurbsCurvePoint(nurbs, u);
        Np= nurbsCurveDerivPoint(nurbs,u);

        xMn1 = x - N(1);
        yMn2 = y - N(2);

        du = sqrt(xMn1^2 + yMn2^2 );
        dpu = -( xMn1*Np(1) + yMn2*Np(2) );

        % Newton-Raphson iteration
        if dpu == 0
            % Small perturbation to avoid division by zero
            u = u + tol;
        else
            u = u - (du^2)/dpu;
        end

        % Control de convergencia
        xMn1 = x - N(1);
        yMn2 = y - N(2);
        du = sqrt(xMn1^2 + yMn2^2 );
        if du<tol     
            Uunique = unique(nurbs.U);
            % if u is close to a knot assign this value
            p = find(abs(Uunique-u)<tolParam);
            if ~isempty(p)
                u = Uunique(p);
            end          
            parameter=u;
            j=nIterMax;
        else
            if j==nIterMax-1
                disp('------> No parameter assigned in PointInvesion!!!!!!!!!');
                parameter = -1;
            end
        end
        j = j+1;
    end
end
end