function sol = postprocessInnerNodesSolution_static(T,elemInfo,faceNodesNum,u,theReferenceElement)

%Number of nodes and elements
innerNodes = theReferenceElement.innerNodes;
nOfInnerNodes = length(innerNodes);
faceNodes = unique(theReferenceElement.faceNodes);
nOfFaceNodes = length(faceNodes);
sol = zeros(length(faceNodesNum),1);
nOfElements = size(T,1);

%Auxiliar vectors
innerVec = 1:nOfInnerNodes;
faceVec = 2:nOfFaceNodes+1;

%Loop in mesh elements
for iElem = 1:nOfElements
    Te = T(iElem,:);
    ifaceNodes = Te(faceNodes);
    inewFaceNodes = faceNodesNum(ifaceNodes);
    faceSol = u(inewFaceNodes);
    
    %Solution at inner nodes
    innerSol = -elemInfo(innerVec,1,iElem) -... %the first minus is because data.mesh.fvolume is stored as positive
        elemInfo(innerVec,faceVec,iElem)*faceSol;
    
    %Complete solution
    sol(Te) = [faceSol ; innerSol];
end
    
    