clear all
close all

addpath(genpath(pwd))

meshFileName = 'mesh1_P1.dcm';
outMeshName = 'mesh1_P4[linear].mat';
refP = 4;
nOfNodesRef = (refP+1)*(refP+2)/2;

meshFileName2Load = GenerateMatFileFromEZ4U_NEFEM(meshFileName);
% meshFileName2Load = {meshFileName};
load(meshFileName2Load{:});
delete(meshFileName2Load{:})

if size(T,2) > 3
    error('Base mesh is not linear!')
end

referenceElement = createReferenceElement(1,nOfNodesRef);
faceNodes = referenceElement.faceNodes(:,2:end-1); %without vertices
coordRef = referenceElement.NodesCoord;
nOfElementNodes = size(coordRef,1);
nOfFaceNodes = size(faceNodes,2);

[intFaces,elemIntfaceInfo] = GetFaces_mod(T);
nOfInteriorFaces = size(intFaces,1);

nOfElements = size(T,1);
nOfNodes = size(X,1);
Xp = zeros(nOfElementNodes*nOfElements,2);
Tp = zeros(nOfElements,nOfElementNodes);
Xp(1:nOfNodes,:) = X;
Tp(:,1:3) = T;

%Mesh
elemPos = [1 3];
facePos = [2 4];
auxiliarCoordLogical = true(nOfElementNodes,1);
auxiliarCoordLogical(1:3) = false;
intFacesMeshed = false(nOfInteriorFaces,1);
ini = nOfNodes + 1;
for elem = 1:nOfElements
    
    coordLogical = auxiliarCoordLogical;
    nOfAlreadyFacesMeshed = 0;
    
    faceInfo = elemIntfaceInfo(elem,:);
    intElemFaces = faceInfo(logical(faceInfo));
    areAlreadyMeshed = intFacesMeshed(intElemFaces,1);
    if any(areAlreadyMeshed)
        facesAlreadyMeshed = intElemFaces(areAlreadyMeshed);
        nOfAlreadyFacesMeshed = length(facesAlreadyMeshed);
        for iface = 1:nOfAlreadyFacesMeshed
            faceMeshed = facesAlreadyMeshed(iface);
            elements = intFaces(faceMeshed,elemPos);
            elementCondition = elements ~= elem;
            elementAlreadyMeshed = elements(elementCondition);
            elementFaceAlreadyMeshed = intFaces(faceMeshed,facePos(elementCondition));
            elementFace = intFaces(faceMeshed,facePos(~elementCondition));
            nodesFaceAlreadyMeshed = faceNodes(elementFaceAlreadyMeshed,:);
            nodesFace = faceNodes(elementFace,:);
            
            Tp(elem,nodesFace) = fliplr(Tp(elementAlreadyMeshed,nodesFaceAlreadyMeshed));
            
            coordLogical(nodesFace) = false;
        end
    end
    
    Te = T(elem,:);
    coordRefMod = coordRef(coordLogical,:);
    nodesElemMod = linearMapping(X(Te,:),coordRefMod);
    
    ind = ini:ini-1 ...
                 + nOfElementNodes ...
                 - nOfAlreadyFacesMeshed*nOfFaceNodes ...
                 - 3; %without vertices
    
    Xp(ind,:) = nodesElemMod;
    Tp(elem,coordLogical) = ind;
    
    intFacesMeshed(intElemFaces) = true;
    if ~isempty(ind)
        ini = ind(end) + 1;
    end

end

Xp(ini:end,:) = [];

%Save
X = Xp;
T = Tp;
elemInfo.nOfNodes = nOfElementNodes;
elemInfo.faceNodes1d = referenceElement.faceNodes1d;
elemInfo.faceNodes = referenceElement.faceNodes;
save(outMeshName,'T','X','Tb_*','elemInfo')
    
    
    
    
    
    
    
    
    
    
    
    



        