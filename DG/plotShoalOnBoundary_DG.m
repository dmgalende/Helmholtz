function plotShoalOnBoundary_DG(X,T,U,referenceElement,infoFaces)
U = reshape(U(:,4,:),numel(T),1);
elemFaceNodes = referenceElement.faceNodes;
nOfElementNodes = size(referenceElement.NodesCoord,1);
nOfCheckFaces = size(infoFaces,1);
faceGNodes = zeros(numel(referenceElement.faceNodes1d),size(infoFaces,1)/2);
refSol = zeros(numel(referenceElement.faceNodes1d),size(infoFaces,1)/2);
j = 1;
for i = 1:nOfCheckFaces
    infoFace = infoFaces(i,:);
    elem = infoFace(1);
    face = infoFace(2);
    nodesElem = T(elem,:); % element nodes in global numbering
    if X(nodesElem(elemFaceNodes(face,1)),2)<1e-5
        faceGNodes(:,j) = nodesElem(elemFaceNodes(infoFace(2),:))'; % face nodes in global numbering
        ind = (elem-1)*nOfElementNodes+1:elem*nOfElementNodes;
        refSol_aux = real(U(ind));
        refSol(:,j) = refSol_aux(elemFaceNodes(infoFace(2),:));
        j = j+1;
    end
end
[ordRefSol,ordFaceGNodes] = concatenateRefSol(refSol,faceGNodes);
refSol = ordRefSol; faceGNodes = ordFaceGNodes;

postprocessRefSolRect(refSol,faceGNodes,X)


function [orderedRefSol,orderedFaceNodes] = concatenateRefSol(refSol,faceNodes)

orderedRefSol = refSol;
orderedFaceNodes = faceNodes;
% find the first face (for open bounaries)
for i = 1:size(faceNodes,2)
    if isempty(find(faceNodes(1,i) == faceNodes(end,:), 1)) %then it's the 1º face
        orderedFaceNodes(:,1) = faceNodes(:,i);
        orderedRefSol(:,1) = refSol(:,i);
        break
    end
end
for i = 1:size(faceNodes,2)-1
    lastFaceNode = orderedFaceNodes(end,i); % find the last node of the i-th face
    suxCol = find(faceNodes(1,:) == lastFaceNode); % find the face to concatenate
    orderedFaceNodes(:,i+1) = faceNodes(:,suxCol);
    orderedRefSol(:,i+1) = refSol(:,suxCol);
end


function postprocessRefSolRect(refSol,faceGNodes,X)
% plot the solution on the reference boundary
% when the boundary is rectilinear
xdata = zeros(size(refSol,2),1);
for i = 1:size(refSol,2)
    cartCoord = X(faceGNodes(:,i),:);
    xdata_rel = coordinate_X(cartCoord(:,1),cartCoord(:,2));
    xdata = xdata(end) + xdata_rel;
    plot(xdata,refSol(:,i),'g')
    hold on
end


function xdata = coordinate_X(x,y)
xdata = zeros(size(x,1),1);
for i = 2:length(xdata)
    xdata(i) = sqrt((x(i)-x(1))^2+(y(i)-y(1))^2);
end
