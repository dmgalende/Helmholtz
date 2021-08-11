function plotSolutionPadaptivity(mesh,u,p_adaptivity,nDegRef)

% Check input
if nargin == 3
    nDegRef = 15;
end

% p refinement variables
elem_p = p_adaptivity.elem_p;
elem_check = p_adaptivity.elem_check;
different_p = p_adaptivity.different_p;
nOfDifferentP = numel(different_p);
allElements = 1:size(mesh.T,1);

% Plotting element (equal spaced points)
nodes = [];
h = 1/nDegRef;
for j = 0:nDegRef
    i = (0:nDegRef-j)';
    aux = j*ones(size(i));
    nodes = [nodes; [i aux]*h];
end
nodes = 2*nodes - 1;
npoints = size(nodes,1);

% Delaunay triangulation of plotting element
elemTriRef = delaunayn(nodes);

% Mesh info
nameNodal = mesh.fieldNames{mesh.indexElemPosCon(2)};
nameCon = mesh.fieldNames{mesh.indexElemPosCon(3)};
X = mesh.(nameNodal);
T = mesh.(nameCon);
nOfElements = size(T,1);

for iP = 1:nOfDifferentP

    if iP==1
        % base elements
        referenceElement = mesh.referenceElement;
        coordRef = referenceElement.NodesCoord;
        nOfNodes = size(coordRef,1);
        nDeg = referenceElement.degree;
        V = Vandermonde_LP(nDeg,coordRef);
        invV = inv(V');
        elements = allElements(~elem_check);
    else
        % refined elements
        referenceElement = p_adaptivity.referenceElements.(['P' num2str(different_p(iP))]);
        coordRef = referenceElement.NodesCoord;
        nOfNodes = size(coordRef,1);
        nDeg = referenceElement.degree;
        V = Vandermonde_LP(nDeg,coordRef);
        invV = inv(V');
        elements = allElements(elem_p==nDeg);
    end
    
    % Compute shape functions at interpolation points
    shapeFunctions = zeros(npoints,nOfNodes);
    for ipoint = 1:npoints
        p = orthopoly2D(nodes(ipoint,:),nDeg);
        shapeFunctions(ipoint,:) = (invV*p)';
    end

    % Loop in elements
    patchHandle = zeros(1,nOfElements);
    for ielem = elements

        % Interpolate solution and position at interpolation points
        Te = T(ielem,:);
        Xe = linearMapping(X(Te(1:3),:),coordRef);
        Xplot = shapeFunctions*Xe;
        ind = 0.5*sum((elem_p(1:ielem-1)+1).*(elem_p(1:ielem-1)+2))+1 : ...
        0.5*sum((elem_p(1:ielem)+1).*(elem_p(1:ielem)+2));
        uplot = shapeFunctions*u(ind);

        % Plot interpolated solution in the element
        hold on
        patchHandle(ielem) = patch('Faces',elemTriRef,'Vertices',Xplot,'FaceVertexCData',uplot,...
            'FaceColor','interp','EdgeAlpha',0);
        hold off
    end
    axis equal
    colorbar
end
