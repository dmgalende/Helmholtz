function data = mesh_preprocess_old(data)

t = cputime;

mesh = data.mesh;
nameCon = mesh.fieldNames{data.mesh.indexElemPosCon(3)};
% Convert Types to use M++
mesh.(nameCon) = int32(mesh.(nameCon));
% Creating interior faces
Tlinear = mesh.(nameCon)(:,1:3);
elements = (1:size(Tlinear,1));
intFaces = GetFaces(Tlinear);
nOfExteriorElements = 0;
% Creating infoFaces structure
nOfBoundaries = numel(data.mesh.boundaryNames);
elementFaceInfo = data.mesh.elementFaceInfo;
elemFaceNodes = data.mesh.referenceElement.faceNodes;
for iboundary = 1:nOfBoundaries
    iname = data.mesh.boundaryNames{iboundary};
    infoFaces.(['exteriorFaces_' iname]) = int32(elementFaceInfo.(iname));
    nOfExteriorElements = nOfExteriorElements + size(mesh.(['Tb_',iname]),1);
end
% Creating zones (int/ext - PML/free-space)
T = mesh.(nameCon);
nOfElements = numel(elements);
nOfInteriorElements = nOfElements - nOfExteriorElements;
interiorElements = 1: nOfInteriorElements;
exteriorElements = setdiff(elements, interiorElements);
elementsPML = [];
interiorFaces_PML = [];
interiorFaces_FreeSpace = intFaces;
candidatesPMLelements = [];
if  any(strcmp(data.PML(5,:),'on'))
    tol = 1e-9;
    % Interior faces are splitted in (1) interior faces in free-space and
    % (2) interior faces in the PML region
    % Different numerical fluxes will be applied in (2) for the Berenger PML
    % Also mark candidate elements to be in the PML region because a source
    % must be defined in this region.
    interiorFaces_FreeSpace = [];
    nOfInteriorFaces = size(intFaces,1);
    for iIntFace = 1:nOfInteriorFaces
        infoFace = intFaces(iIntFace,:);
        iElem = infoFace(1);
        nodesElem = T(iElem,:); % element nodes in global numbering
        nodesFace = nodesElem(elemFaceNodes(infoFace(2),:)); % face nodes in global numbering
        absValueFace = data.PMLabsorptionValue(nodesFace,:); % absorption value in the face nodes
        absValueFace(absValueFace<tol) = 0; 
        if all(absValueFace(:,1)~=0) || all(absValueFace(:,2)~=0)
            interiorFaces_PML = [interiorFaces_PML; infoFace];
            candidatesPMLelements = [candidatesPMLelements, infoFace([1,3])];
        else
            interiorFaces_FreeSpace = [interiorFaces_FreeSpace; infoFace];
        end
    end
    candidatesPMLelements = unique([candidatesPMLelements exteriorElements]);
    % Elements on the PML region are identified. To reduce the computational
    % cost of this loop candidate elements are obtained in the previos loop.
    % Candidate elements are elements having a face in the PML region.
    for iElem = candidatesPMLelements
        nodes = T(iElem,:);
        absValueElem = data.PMLabsorptionValue(nodes,:);
        absValueElem(absValueElem<tol) = 0; 
        if all(absValueElem(:,1)~=0) || all(absValueElem(:,2)~=0)
            elementsPML = [elementsPML, iElem];
        end
    end
    elementsPML = unique(elementsPML);
end
interiorElements_FreeSpace = int32(setdiff(interiorElements, elementsPML));
interiorElements_PML = int32(setdiff(interiorElements, interiorElements_FreeSpace));
exteriorElements_FreeSpace = int32(setdiff(exteriorElements, elementsPML));
exteriorElements_PML = int32(setdiff(exteriorElements, exteriorElements_FreeSpace));
% Store information--------------------------------------------------------
infoFaces.interiorFaces_FreeSpace = int32(interiorFaces_FreeSpace);
infoFaces.interiorFaces_PML = int32(interiorFaces_PML);
mesh.interiorElements_FreeSpace =interiorElements_FreeSpace;
mesh.interiorElements_PML = interiorElements_PML;
mesh.exteriorElements_FreeSpace = exteriorElements_FreeSpace;
mesh.exteriorElements_PML = exteriorElements_PML;

data.infoFaces = infoFaces;
data.mesh = mesh;
data.cputime.preprocess_DG = cputime - t;
