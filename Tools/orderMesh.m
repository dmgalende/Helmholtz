function p = orderMesh(mesh, boundary, byelements)

if nargin == 2, byelements = false; end 

% Some variables
nodes = unique(mesh.T);
info = mesh.elementFaceInfo.(boundary);
elems_boundary = info(:, 1);

% Boundary nodes
if byelements                               % use nodes from elements having edges on the boundary
    nodes_boundary = unique(mesh.T(elems_boundary, :));
else                                        % use nodes strictly on the boundary
    face_nodes = mesh.refelem.faceNodes;
    faces_boundary = info(:, 2);
    nElems = length(elems_boundary);
    nodes_boundary = zeros(nElems, size(face_nodes, 2));
    for i = 1:nElems
        nodes_boundary(i, :) = mesh.T(elems_boundary(i), face_nodes(faces_boundary(i), :));    
    end
    nodes_boundary = unique(nodes_boundary);
end

% Interior nodes
nodes_rest = setdiff(nodes, nodes_boundary);
ref = length(nodes_rest);

% Mapping to order DOFs in boundary elements
p = nodes;
p(nodes_boundary) = (1:length(nodes_boundary)) + ref;

% Mapping to order DOFs not in boundary elements
p(nodes_rest) = 1:ref;
