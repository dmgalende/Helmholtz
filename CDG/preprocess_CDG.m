function data = preprocess_CDG(data)
t = cputime;
% Creating interior faces
nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};
Tlinear = data.mesh.(nameCon)(:,1:3);    
intFaces = GetFaces(Tlinear);
infoFaces.interiorFaces = intFaces;

% Creating infoFaces structure
nOfBoundaries = numel(data.mesh.boundaryNames);
elementFaceInfo = data.mesh.elementFaceInfo;
for iboundary = 1:nOfBoundaries
    iname = data.mesh.boundaryNames{iboundary};
    infoFaces.(['exteriorFaces_' iname]) = elementFaceInfo.(iname);
end
data.infoFaces = infoFaces;
data.cputime.infoFaces = cputime - t;