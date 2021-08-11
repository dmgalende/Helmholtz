function [gaussPoints,weightsPoints] = nefemQuad2DElementLocalCoordinates...
                            (vertCoord,aNurbs,u1,u2,nIPNurbs,nIPinterior)
%
% [gaussPoints,weightsPoints] = computeQuadNurbsElem...
%               (vertCoord,aNurbs,u1,u2,nIPNurbs,nIPinterior);
%
% Quadrature definition for a nurbs element
% OBS: The changes of definition of the nurbs boundary have 
%      kept in mind itself
%
% Input:
% vertCoord:     nodal coordinates of the element
% aNurbs:        description of the nurbs boundary 
% u1,u2:         parameters of the nurbs curve (trimmed nurbs)
% nIPNurbs, nIPinterior:number of integration points in u (v) direction
%
% Output:
% gaussPoints, weightsPoints: integration points and weights
%

[zU,wU] = gaussLegendre(nIPNurbs,-1,1);
[zV,wV] = gaussLegendre(nIPinterior,-1,1);

% Definition of the quadrature over [-1,1]^2
z = zeros(nIPNurbs*nIPinterior, 2);
w = zeros(1,nIPNurbs*nIPinterior);
k = 1;
for i = 1:nIPNurbs
    for j = 1:nIPinterior
        z(k,:) = [zU(i),zV(j)];
        w(k) = wU(i)*wV(j);
        k = k+1;
    end
end
clear zU zV wU wV;

% Compound quadrature 
Mu = 1;
Mv = 1;

zComp = [];
wComp = [];

hu = 2/Mu;
hv = 2/Mv;
for i=1:Mu
    xIni = -1 + (i-1)*hu;
    xFin = -1 + i*hu;
    for j=1:Mv
        yIni = -1 + (j-1)*hv;
        yFin = -1 + j*hv;

        x = ((z(:,1)+1)/2)*(xFin-xIni) + xIni;
        y = ((z(:,2)+1)/2)*(yFin-yIni) + yIni;

        pes = w*(xFin-xIni)*(yFin-yIni)/4;
        zComp = [zComp; [x y]];
        wComp = [wComp, pes];
    end
end
clear z w;

% Quadrature over [0,1]^2
z01 = (zComp + 1)/2;
w01 = wComp/4;
clear zComp wComp;

% Quadrature over  [u1 u2]x[0 1] (trimmed Nurbs) taking into account 
% the changes of nurbs definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Searching for changes of definition in the NURBS boundary 
% (for the current element)

knotsNurbs = aNurbs.U;
knotsChanges = unique(knotsNurbs(find(knotsNurbs>=min(u1,u2)...
                                    & knotsNurbs<=max(u1,u2))));
knotsQuadrature = unique([min(u1,u2) knotsChanges max(u1,u2)]);
nKnotsQuadrature = length(knotsQuadrature);

z = [];
w = [];
for iKnotsQuadrature = 1:nKnotsQuadrature-1
    zAux1 = z01(:,1)*knotsQuadrature(iKnotsQuadrature+1) ...
                 + (1-z01(:,1))*knotsQuadrature(iKnotsQuadrature);
    zAux2 = z01(:,2);
    wAux = w01*abs(knotsQuadrature(iKnotsQuadrature+1)...
        -knotsQuadrature(iKnotsQuadrature));
    z = [z; zAux1 zAux2];
    w = [w, wAux];
end

% Quadrature over the triangular element with an edge at the aNurbs boundary
nOfGauss = size(z,1);
gaussPoints = zeros(nOfGauss,2);
weightsPoints = zeros(1,nOfGauss);

% Vertex used to parametrize the element
vertexParam = [0 1];

for iGauss = 1:nOfGauss
    u = z(iGauss,1);
    v = z(iGauss,2);
    
    % Over aNurbs curve (local coordinates)
    pt = nurbsCurvePoint(aNurbs,u);
    
    % Inverse of the isoparametric transformation
    [phi01,JLinear] = inverseLinearMapping(vertCoord,pt(1:2));
    phi = (phi01+1)/2;
    
    % Element parametrization
    gaussPoints(iGauss,:) = (1-v)*phi  + v*vertexParam;
    
    % Change of the integration weights
    Fu = nurbsCurveDerivPoint(aNurbs,u);    
    fip = JLinear\[Fu(1);Fu(2)];    
    JParametrization = abs( (1-v)*fip(1)*(1-phi(2)) + phi(1)*(1-v)*fip(2) );    
    weightsPoints(iGauss) = w(iGauss)*JParametrization;    
end

gaussPoints = 2*gaussPoints-1;
weightsPoints = weightsPoints*4;


