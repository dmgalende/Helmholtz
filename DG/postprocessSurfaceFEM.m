function postprocessSurfaceFEM(theMesh, theReferenceElement, U, iComp)

nDegPlot = 20;
optBoundary = 0;

X = theMesh.X;
T = theMesh.T;
nOfElements = size(T,1);
 
if nargin<4
    iComp = 1;
end

nodesRefTri = equallySpacedNodesRefTri(nDegPlot);
connecNodesTri = delaunay(nodesRefTri(:,1), nodesRefTri(:,2));

shapeFunctionsFEM = plotInfo_FEM(theReferenceElement,nDegPlot);

% Loop on elements
for iElem = 1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    
    Ue = U(:,iComp,iElem);
    interpNodesElem = shapeFunctionsFEM'*Xe;
    interpSolElem = shapeFunctionsFEM'*Ue;

    trisurf(connecNodesTri, interpNodesElem(:,1), interpNodesElem(:,2), interpNodesElem(:,2)*0,...
        interpSolElem,'EdgeColor','none','FaceColor','interp','FaceLighting','phong');
    hold on
    if optBoundary
        edgeNodesOrdered = postprocessComputeEdgeNodesOrdered(nDegPlot);
        plot3(interpNodesElem(edgeNodesOrdered,1), interpNodesElem(edgeNodesOrdered,2), ...
            interpSolElem(edgeNodesOrdered,1)*0,'k-','LineWidth',1)
    end
end
view(0,90)
axis off
axis equal