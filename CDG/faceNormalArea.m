% function [n,A] = faceNormalArea(infoFace,theReferenceElement,theMesh)
function [n,A] = faceNormalArea(infoFace,theReferenceElement,X,T)
%
% [n,A] = faceNormalArea(infoFace,theReferenceElement,theMesh)
%
% Function to compute the unitari normal vector (outward respect
% to the first element) and the area of the face
%
% Input:
% infoFace:             information of the current face
% theReferenceElement:  struct containing the reference element
% theMesh:              struct containing the information of the mesh
%
% Output:
% n: unit outward normal
% A: area of the face
%

nElem = infoFace(1);
nFace = infoFace(2);

faceNodes = theReferenceElement.faceNodes;
faceVertices = faceNodes(:,[1 end]);
nodes = T(nElem,faceVertices(nFace,:));
coord = X(nodes,:);

nsd = size(coord,2);
switch nsd
    case 3
        u = coord(2,:)-coord(1,:);   v = coord(3,:)-coord(1,:);
        w = cross(u,v);
        normw = norm(w);
        A = normw/2;
        n = w/normw;
    case 2
        u = coord(2,:)-coord(1,:);
        A = norm(u);
        n = [u(2),-u(1)]/A;
    otherwise
        error('wrong nsd in coordVert')
end