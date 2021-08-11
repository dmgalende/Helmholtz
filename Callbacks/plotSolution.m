function varargout = plotSolution(X,T,u,referenceElement,flag)

% [patchHandle,tri] = plotSolution(X,T,u,referenceElement) plots a nodal 
% scalar field using delaunay triangulation.
% 
% Input:
%   X: nodal coordinates
%   T: FEM connectivity matrix.
%   u: nodal values
%   referenceElement: the reference element structure. 
%
%   NOTE: if referenceElement argument is not given T is assumed as linear, 
%         so triangulation wont be done.
%
% Output:
%   patchHandle (optional): handle to the created patch object
%   tri (optional): linear conectivity matrix for plotting

%Check degree of the element
nOfNodesPerElement = size(T,2);

%Delaunay's solution mesh
if exist('referenceElement','var') && nOfNodesPerElement > 4
    
    % Creating solution mesh conectivity (tri)
    coordRef = referenceElement.NodesCoord;
    elemTriRef = delaunayn(coordRef);
    nOfElemTriRef = size(elemTriRef,1);
    nOfElements = size(T,1);
    tri = zeros(nOfElemTriRef*nOfElements,3);
    indexElem = 0;
    for ielem = 1:nOfElements
        Te = T(ielem,:);
        for ielemRef = 1:nOfElemTriRef
            indexElemRef = indexElem + ielemRef;
            tri(indexElemRef,:) = Te(elemTriRef(ielemRef,:));
        end
        indexElem = indexElem + nOfElemTriRef;
    end
    
else 
    tri = T;
end

% Return tri if flag argument is given (no plot)
if exist('flag','var') && flag
    varargout = {tri};
    return
end

% Plot
patchHandle = patch('Faces',tri,'Vertices',X,'FaceVertexCData',u,...
    'FaceColor','interp','EdgeAlpha',0);
axis equal
colorbarHandle = colorbar('location','East');

%Output variable
if ~nargout
    varargout = {[]};
elseif nargout == 1
    varargout = {[patchHandle colorbarHandle]};
elseif nargout == 2
    varargout = {[patchHandle colorbarHandle] tri};
end
