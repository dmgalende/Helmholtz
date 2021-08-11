function MatRefElem = computeReferenceElementMatrices(referenceElement)
%
% MatRefElem = computeReferenceElementMatrices(referenceElement)
%
% Function to compute all the elemental matrices for the reference element
%
% Input: 
% referenceElement: struct containing the information of the reference element 
%
% Output:
% MatRefElem: struct containing the elemental matrices for the reference
%             element
%

 n = referenceElement.degree;
 coord = referenceElement.NodesCoord;
 
%Vandermonde matrix
 V = Vandermonde_LP(n,coord);

%Inverse of the mass matrix
 invM = V*V';
 
%Inverse of the Vandermonde matrix
 invV = inv(V);
 
 nsd = size(coord,2);

%  MatRefElem.nOfSpatialDimensions = nsd;
%  MatRefElem.interpolationDegree = n;
 MatRefElem.invM = invM;
%  MatRefElem.M = invV'*invV;
 MatRefElem.invV = invV;

 %Convection matrices
 switch nsd
     case 1,
         Cxi = computeReferenceElementConvectionMatrices(n,invV);
         MatRefElem.Cxi = Cxi;
     case 2,
         [Cxi,Ceta] = computeReferenceElementConvectionMatrices(n,invV);
         MatRefElem.Cxi = Cxi;        MatRefElem.Ceta = Ceta;
     case 3,
         [Cxi,Ceta,Czeta] = computeReferenceElementConvectionMatrices(n,invV);
         MatRefElem.Cxi = Cxi;      MatRefElem.Ceta = Ceta;      MatRefElem.Czet = Czeta;
     otherwise
         error('computeReferenceElementMatrices can only be used in 1D, 2D or 3D')
 end