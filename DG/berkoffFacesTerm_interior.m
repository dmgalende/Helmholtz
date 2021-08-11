function F = berkoffFacesTerm_Interior(X,T,U,theReferenceElement,infoFaces,faceMassReference,ccg)

F = zeros(size(U));
nOfNodes = size(U,1);
nOfComponents = size(U,2);
elemFaceNodes = theReferenceElement.faceNodes;

%% Interior faces in free-space
nOfInteriorFaces_FreeSpace = size(infoFaces.interiorFaces_FreeSpace,1);
for i = 1:nOfInteriorFaces_FreeSpace
    infoFace = infoFaces.interiorFaces_FreeSpace(i,:);    
    iElemL = infoFace(1);
    iElemR = infoFace(3);
    [n,A] = faceNormalArea(infoFace,theReferenceElement,X,T); 
    % bottom informations
    nodesElem = T(iElemL,:); % element nodes in global numbering
    faceLocalNodes_L = elemFaceNodes(infoFace(2),:); % face nodes in local numbering Left Elem
    faceLocalNodes_R = elemFaceNodes(infoFace(4),end:-1:1); % face nodes in local numbering Right Elem
    faceNodes = nodesElem(faceLocalNodes_L); % face nodes in global numbering Left Elem
    ccg_face = ccg(faceNodes);
    % Left and right states
    uL = U(faceLocalNodes_L,:,iElemL);
    uR = U(faceLocalNodes_R,:,iElemR);
    % Numerical minus normal flux-------------------
    [FluxUL,FluxUR] = berkoffFluxU(uL, uR, n, ccg_face);    
    % Left element----------------------------------
    Fe = zeros(nOfNodes,nOfComponents);
    Ff = faceMassReference*FluxUL*A/2;
    Fe(faceLocalNodes_L,:) = Ff;
    F(:,:,iElemL) = F(:,:,iElemL) + Fe;    
    % Right element----------------------------------
    Fe = zeros(nOfNodes,nOfComponents);
    Ff = faceMassReference*FluxUR*A/2;
    Fe(faceLocalNodes_R,:) = Ff;
    F(:,:,iElemR) = F(:,:,iElemR) + Fe;
end

%% Interior faces in PML
nOfInteriorFaces_PML = size(infoFaces.interiorFaces_PML,1);
for i = 1:nOfInteriorFaces_PML
    infoFace = infoFaces.interiorFaces_PML(i,:);    
    iElemL = infoFace(1);
    iElemR = infoFace(3);
    [n,A] = faceNormalArea(infoFace,theReferenceElement,X,T); 
    % bottom informations
    nodesElem = T(iElemL,:); % element nodes in global numbering
    faceLocalNodes_L = elemFaceNodes(infoFace(2),:); % face nodes in local numbering Left Elem
    faceLocalNodes_R = elemFaceNodes(infoFace(4),end:-1:1); % face nodes in local numbering Right Elem
    faceNodes = nodesElem(faceLocalNodes_L); % face nodes in global numbering
    ccg_face = ccg(faceNodes);
    % Left and right states
    uL = U(faceLocalNodes_L,:,iElemL);
    uR = U(faceLocalNodes_R,:,iElemR);
    % Numerical minus normal flux-------------------
    [FluxUL,FluxUR] = berkoffPMLFluxU(uL, uR, n, ccg_face);    
    % Left element----------------------------------
    Fe = zeros(nOfNodes,nOfComponents);
    Ff = faceMassReference*FluxUL*A/2;
    Fe(faceLocalNodes_L,:) = Ff;
    F(:,:,iElemL) = F(:,:,iElemL) + Fe;    
    % Right element----------------------------------
    Fe = zeros(nOfNodes,nOfComponents);
    Ff = faceMassReference*FluxUR*A/2;
    Fe(faceLocalNodes_R,:) = Ff;
    F(:,:,iElemR) = F(:,:,iElemR) + Fe;
end
