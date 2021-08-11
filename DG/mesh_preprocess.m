function data = mesh_preprocess(data)

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
PMLnodes = [];
for iboundary = 1:nOfBoundaries
    iname = data.mesh.boundaryNames{iboundary};
    infoFaces.(['exteriorFaces_' iname]) = int32(elementFaceInfo.(iname));
    nOfExteriorElements = nOfExteriorElements + size(mesh.(['Tb_',iname]),1);
    if strcmp(data.PML{5,iboundary},'on')
        PMLnodes = unique([PMLnodes data.PML{4,iboundary}{1} data.PML{4,iboundary}{2}]);
    end
end
% Creating zones (int/ext - PML/free-space)
T = mesh.(nameCon);
nOfElements = numel(elements);
nOfInteriorElements = nOfElements - nOfExteriorElements;
interiorElements = 1: nOfInteriorElements;
exteriorElements = setdiff(elements, interiorElements);
elementsPML = zeros(nOfElements,1);
if  any(strcmp(data.PML(5,:),'on'))
    T_lin = T(:,1:3);
    aux = ismember(T_lin,PMLnodes);
    elementsPML = elements(all(aux'));
    ind1 = ismember(intFaces(:,1),elementsPML);
    ind2 = ismember(intFaces(:,3),elementsPML);
    interiorFaces_FreeSpace = intFaces(~(ind1|ind2),:);
    interiorFaces_PML = intFaces(ind1|ind2,:);
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
