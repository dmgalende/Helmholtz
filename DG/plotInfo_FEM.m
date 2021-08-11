function shapeFunctionsFEM = plotInfo_FEM(theReferenceElement,nDegPlot)
%
%
% Triangulation of a reference element
nodesRefTri = equallySpacedNodesRefTri(nDegPlot);

nDeg = theReferenceElement.degree;
coordRef = theReferenceElement.NodesCoord;
nOfElementNodes = size(coordRef,1);
V = Vandermonde_LP(nDeg,coordRef);
invVt = inv(V');

% Precomputations for FEM
nNodesRefTri = size(nodesRefTri,1);
shapeFunctionsFEM = zeros(nOfElementNodes,nNodesRefTri);
for iNode = 1:nNodesRefTri
    pRST = orthopoly2D(nodesRefTri(iNode,:),nDeg);
    % Shape function (for all element nodes evaluated on the iPoint)
    shapeFunctionsFEM(:,iNode) = (invVt*pRST)';
end
