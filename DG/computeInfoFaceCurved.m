function quadCurvedFace = computeInfoFaceCurved(theReferenceElement,coordFace)
%
% quadCurvedFace = computeInfoFaceCurved_FEM(theReferenceElement,coordFace)
%
% Store for every curved face
% quadrature on the curved face
% the outward unit normal as gauss points
% shape functions at gauss points
%

gaussWeights = theReferenceElement.IPweights1d;
nOfGauss = length(gaussWeights);

% Face information
quadCurvedFace.weights = zeros(nOfGauss,1);
quadCurvedFace.extNormal = zeros(nOfGauss,2);
quadCurvedFace.gaussPoints = zeros(nOfGauss,2);

N = theReferenceElement.N1d;
Nxi  = theReferenceElement.N1dxi;

x = coordFace(:,1);
y = coordFace(:,2);

%Loop on Gauss points
for iGauss=1:nOfGauss
    n = [Nxi(iGauss,:)*y, -Nxi(iGauss,:)*x];
    normN = norm(n);
    n = n/normN;

    quadCurvedFace.gaussPoints(iGauss,:) = [N(iGauss,:)*x   N(iGauss,:)*y];
    quadCurvedFace.extNormal(iGauss,:) = n;
    quadCurvedFace.weights(iGauss) = normN*gaussWeights(iGauss);
end          
