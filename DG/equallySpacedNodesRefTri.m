function nodes = equallySpacedNodesRefTri(nDeg)
%
%  nodes = equallySpacedNodesRefTri(nDeg)
%  
% Input:
% nDeg: degree of the finite elements
%
% Output:
% nodes: Equally-space distribution of nodes in the 
%        reference triangle [-1,-1; 1,-1; -1,1]
%

h = 1/nDeg;

% Equally-space nodes 
nodes = [];
for j=0:nDeg    
    i = [0:nDeg-j]';
    aux = j*ones(size(i));
    nodes = [nodes; i*h aux*h];
end

% In the reference triangle [-1,-1; 1,-1; -1,1]
nodes = 2*nodes-1;