function data = matrix_preprocess(data,handles)

t = cputime;
nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};
X = data.mesh.(nameNodal);
T = data.mesh.(nameCon);
referenceElement = data.mesh.referenceElement;
infoFaces = data.infoFaces;
exteriorElements = union(data.mesh.exteriorElements_FreeSpace,data.mesh.exteriorElements_PML);
c = data.bottom.c;
cg = data.bottom.cg;
incidentWave = data.ip;
boundaryNames = data.mesh.boundaryNames;
elementFaceInfo = data.mesh.elementFaceInfo;

referenceElement.faceNodes = int32(referenceElement.faceNodes);
referenceElement.faceNodes1d = int32(referenceElement.faceNodes1d);
referenceElement.innerNodes = int32(referenceElement.innerNodes);

setOutput({'Compute elemental matrices for the reference element'},handles);
elementalMatricesReference = computeReferenceElementMatrices(referenceElement);

setOutput({'Compute elemental mass matrix for the reference face...'},handles);
faceMassReference = computeReferenceElementMassMatrix(referenceElement.degree,referenceElement.NodesCoord1d);

setOutput({'Compute elemental matrices for exterior elements'},handles)
elementalMatricesExterior = computeAllElementalMatrices(X,T(exteriorElements,:),referenceElement);

setOutput({'Compute elemental vectors for all elements'},handles);
elementalVectors = computeAllElementalVectors(X,T,referenceElement,c,cg,incidentWave);

setOutput({'Compute information for all curved faces...'},handles)
ExteriorFaces = computeAllFacesInfo(X,T,referenceElement,infoFaces,boundaryNames,elementFaceInfo);

elementalMatricesInfo.elementalMatricesReference =  elementalMatricesReference;
elementalMatricesInfo.faceMassReference = faceMassReference;
elementalMatricesInfo.elementalMatricesExterior = elementalMatricesExterior;
elementalMatricesInfo.ExteriorFaces = ExteriorFaces;
elementalMatricesInfo.elementalVectors = elementalVectors;

data.mesh.referenceElement = referenceElement;
data.BC.values = int32(data.BC.values);
data.mesh.indexElemPosCon = int32(data.mesh.indexElemPosCon);
data.elementalMatricesInfo = elementalMatricesInfo;
data.cputime.matrix_preprocess = cputime - t;