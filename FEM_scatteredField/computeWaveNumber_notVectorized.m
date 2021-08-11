function k = computeWaveNumber(omega,bottom)

%Number of nodes
nOfNodes = size(bottom,1);

%Memory allocation
k = zeros(nOfNodes,1);

%Loop in nodes
for inode = 1:nOfNodes
    h = bottom(inode);
    k(inode) = evaluateDispersionRelation(omega,h);
end
k = abs(k);

%_______________________________________
function k = evaluateDispersionRelation(arg,h)

% The function evaluateDispersionRelation(arg,h,varargin) solves the 
% non linear dispersion equation in terms of the wave number k:
%
%                               w^2 = k*g*tanh(k*h)
%
% Input:
%   arg: could be a vector [omega k0] where omega is the angular frequency  
%        and k0 is the initial aproximation of k. If this last value is not
%        given, a default one k0 = omega/sqrt(g*h) will be considered.
%   h: vector with the bottom values for each node.
% Output:
%   k: wave number result of solving dispersion equation

%Initial guess
omega = arg(1);
if ~isscalar(arg), k0 = arg(2); else k0 = omega/sqrt(9.81*h); end

%Problem definition and solver
dispersion = @(k)omega^2-k*9.81*tanh(k*h);
options = optimset('Display','off');
[k,dispersionValue,flag] = fzero(dispersion,k0,options);

%Display warnings
switch flag
    case -3
        warning('NaN or Inf detected in evaluateDispersionRelation') %#ok<WNTAG>
    case -4
        warning('Complex dispersion value detected in evaluateDispersionRelation') %#ok<WNTAG>
    case -5
        warning('evaluateDispersionRelation may have converged to a singular point') %#ok<WNTAG>
end