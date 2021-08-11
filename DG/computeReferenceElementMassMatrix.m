function ElemM = computeReferenceElementMassMatrix(nDeg,coord)
%
% ElemM = computeReferenceElementMassMatrix(nDeg,coord)
%
% Function to compute the elemental mass matrix for the reference element
%
% Input:
% nDeg:   degree of polynomials
% coord:  nodal coordinates (size is nOfNodes X nsd)
% 
% Output:
% ElemM:      Mass matrix 
%             size is nOfNodes X nOfNodes
%

% Vandermonde matrix
invV = inv(Vandermonde_LP(nDeg,coord));

% Mass matrix
ElemM = invV'*invV;