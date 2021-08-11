function [gauss,weight,extNormGauss] = nefemQuad2DFaceParametricSpace(aNurbs, u1, u2, nIP)
    
%
% [gauss,weight,extNormGauss] = nefemQuad2DFaceParametricSpace(aNurbs, u1, u2, nIP)
%
% Quadrature definition for a nurbs boundary integral ("nurbs face")
% in cartesian coordinates
% OBS: The changes of definition of the nurbs have kept in mind itself
%
% Input:
% aNurbs:     description of the nurbs boundary 
% u1,u2:      parameters of the aNurbs curve (trimmed nurbs)
% nIP:        number of integration points
%
% Output:
% gauss, weight: integration points and weights (physical coordinates)
% extNormGauss:  outward unit normal
%

u1Old = u1;
u2Old = u2;

aux = min(u1, u2);
u2 = max(u1, u2);
u1 = aux;

% Search for breakpoints in the interval uInt______________________
uInt = u1;
breakpoints = unique(aNurbs.U);
nBreakpoints = length(breakpoints);
for iBreakpoints = 1:nBreakpoints
    iBreak = breakpoints(iBreakpoints);    
    if iBreak>u1 && iBreak<u2
        uInt = [uInt, iBreak];
    end
end
uInt = [uInt, u2];

if u1Old > u2Old
    uInt = fliplr(uInt);
end
    
% Number of subintervals__________________________________________
nInt = size(uInt,2) - 1;

%  Quadrature definition in [-1,1]
[pospg,pespg]=gaussLegendre(nIP,-1,1);

gauss = [];
weight = [];
extNormGauss = [];

% Loop over subintervals (changes of definition!)_________________
% There are no breakpoints in the subintervals.___________________
for i = 1:nInt
    uIni = uInt(i);
    uFin = uInt(i+1);
    Lu = uFin-uIni;

    % Quadrature definition in the parametric space
    pospgu = (pospg+1)*(Lu/2) + uIni;
    pespgu = pespg*(abs(Lu)/2);
    gauss = [gauss, pospgu];
    weight = [weight, pespgu];
end