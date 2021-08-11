function varargout = plotMesh(X,T,faceNodes,option,nodesNum)

% plotMesh(X,T,faceNodes,option,nodesNum) plots the mesh defined by X and T
% 
% Input:
%   X: nodal coordinates
%   T: connectivities (elements)
%   faceNodes: nOfFaces x nOfNodesPerFace matrix. Indicates the order of
%              the face nodes in a reference element. The numbering of
%              the faces has to be given in a clockwise sense, for instance:
%
%                       QUA elements                 TRI elements     
%
%                             3
%                         *-------*                       *                
%                         |       |                      / \
%                       4 |       | 2                 3 /   \ 2            
%                         |       |                    /     \               
%                         *-------*                   *-------* 
%                             1                           1
%
%              For a given face, the column index of the matrix indicates 
%              the global position of the node. This global numbering has
%              to ascend in a clockwise sense too.
%
%   option (optional): type 'plotNodes' to see the nodes' position on the 
%                      ploted mesh, or type 'plotNodesNum' to see their global 
%                      number.
%   nodesNum (necesary if option = 'plotNodesNum'): type 'all' to plot the
%                                                   global postion of all
%                                                   nodes, or enter a list
%                                                   array with the selected
%                                                   nodes.
%
% Output:
%   patchHandle (optional): handle to the created patch object


%Ordering the face nodes in a row vector without connectivity between them
[nOfFaces,nOfNodesPerFace] = size(faceNodes);
oFaceNodes = zeros(1,nOfFaces*(nOfNodesPerFace-1));
np = nOfNodesPerFace - 1;
aux = 1 - np;
aux2 = 0;
for iface = 1:nOfFaces
    aux = aux + np;
    aux2 = aux2 + np;
    oFaceNodes(aux:aux2) = faceNodes(iface,1:np);
end

%Conectivity for the faces
patchFaces = T(:,oFaceNodes);

%Plot mesh
patchHandle = patch('Faces',patchFaces,'Vertices',X,'FaceColor',[1 1 1],'EdgeAlpha',1);
axis equal

%Optional plots
nodes = unique(T);
if nargin > 3
    hold on
    if strcmpi(option,'plotNodes')
        plot(X(nodes,1),X(nodes,2),'o','markerSize',3,'markerFaceColor','b')
    elseif (nargin == 5) && strcmpi(option,'plotNodesNum')
        if strcmpi(nodesNum,'all')
            list = nodes;
            fontSize = 10;
        elseif ~isnumeric(nodesNum)
            error('wrong list of nodes')
        else
            list = nodesNum;
            fontSize = 15;
            plot(X(list,1),X(list,2),'o','markerSize',3,'markerFaceColor','b')
        end
        for inode = list
            text(X(inode,1),X(inode,2),int2str(inode),'FontSize',fontSize,...
                'Color',[1 0 0])
        end
    elseif (nargin == 4) && strcmpi(option,'plotElemNum')
        for elem = 1:size(T,1);
            xp = mean(X(T(elem,:),1));
            yp = mean(X(T(elem,:),2));
            text(xp,yp,num2str(elem));
        end
    else
        error('wrong optional argument. Check help to fix the error')
    end
    hold off
end

%Output variable
if ~nargout
    varargout = {[]};
else
    varargout = {patchHandle};
end