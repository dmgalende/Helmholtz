function ExteriorFaces = computeAllFacesInfo(X,T,referenceElement,infoFaces,boundaryNames,elementFaceInfo)

nOfBoundaries = numel(boundaryNames);

for iname = 1:nOfBoundaries
    boundary =  boundaryNames{iname};
    infoFacesAux = eval(['infoFaces.exteriorFaces_',boundary]);
    nExteriorFace = size(infoFacesAux,1);
    for i = 1:nExteriorFace
        iElem = infoFacesAux(i,1);
        iFace = infoFacesAux(i,2);
        nodesElem = T(iElem,:); % element nodes in global numbering
        locFaceNodes = referenceElement.faceNodes(iFace,:); % face nodes local numbering
        globFaceNodes = nodesElem(locFaceNodes); % face nodes global numbering
        coordFace = X(globFaceNodes,:);
        ExteriorFacesAux(i) =  computeInfoFaceCurved(referenceElement,coordFace);
    end
    if exist('ExteriorFacesAux','var')
        eval(['ExteriorFaces.',boundary,'= ExteriorFacesAux;']);
    else
        eval(['ExteriorFaces.',boundary,' = [];']);
    end
    clear ExteriorFacesAux;
end
