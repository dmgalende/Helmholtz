function varargout = plotSolution_NEFEM(mesh,u,varargin)

%
% [patchHandle,tri] = plotSolution_NEFEM(mesh,u,T_linear,npoints) plots a 
% nodal scalar field for a NEFEM input mesh using delaunay triangulation.
% 
% Input:
%   mesh: structure containing the X, T, Tb_, matrices and their 
%         corresponding indexes and name in the structure. Also the nurbs, 
%         the trimmedInfo and the reference element have to be inclouded.
%   u: nodal values
%   T_linear (optional): linear connectivity matrix for plotting interior 
%                        elements. If it is given the triangulation for
%                        these elements wont be done.
%   nDegRef (optional): reference degree for the plot element. 
%                       By default is 3.
%
%   NOTE: this function calls plotSolution.m for plotting interior elements.    
%
% Output:
%   patchHandle (optional): handle (vector of handles) to the patch objects
%   tri (optional): linear conectivity matrix for plotting interior
%                   elements (same as T_linear, if it is given).
%

% Check optional input arguments
if length(varargin) == 1
    arg = varargin{:};
    if isscalar(arg)
        nDegRef = arg;
        checkT_linear = false;
    else
        nDegRef = 20;
        tri = arg;
        checkT_linear = true;
    end
elseif length(varargin) == 2
    arg1 = varargin{1};
    arg2 = varargin{2};
    checkT_linear = true;
    if isscalar(arg1)
        nDegRef = arg1;
        tri = arg2;
    else
        nDegRef = arg2;
        tri = arg1;
    end
elseif isempty(varargin)
    checkT_linear = false;
    nDegRef = 3;
else
    error('error with input arguments')
end

% Plot interior elements
if checkT_linear
    h = plotSolution(mesh.X,tri,u);
else
    h = plotSolutionInterpElem(mesh,u,'FEM',nDegRef,'T');
end
delete(h(end))
patchHandle = h(1:end-1);

%Info
referenceElement = mesh.referenceElement;
nDeg = referenceElement.degree;
coordRef = referenceElement.NodesCoord;
nOfNodes = size(coordRef,1);

% Nodes for plot
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
%elemTriRef = delaunayn(nodes,{'QJ','QbB','Qc','Qx','Pp'});

% Plot boundary elements
nurbs = mesh.nurbs;
nOfBoundaries = numel(mesh.boundaryNames);
for cont = 1:nOfBoundaries
    item = mesh.boundaryIndex(cont);
    iboundary = mesh.fieldNames{item};
    iboundaryTrimInfo = mesh.trimmedInfo.(iboundary);
    iboundaryData = mesh.(iboundary);
    nOfElements = size(iboundaryData,1);
    elemPatchHandle = zeros(1,nOfElements);

    % Loop in boundary elements
    for ielem = 1:nOfElements

        % Element connectivity
        Te = iboundaryData(ielem,:);

        % NEFEM info
        u1 = iboundaryTrimInfo(ielem).trim(1);
        u2 = iboundaryTrimInfo(ielem).trim(2);
        idNurbs = iboundaryTrimInfo(ielem).idNurbs;
        aNurbs = nurbs(idNurbs);
        vertCoord = mesh.X(iboundaryData(ielem,1:3),:);

        % Adapted reference element nodes to nurbs size
        nodesAdapted = nefemInterp2DAdaptedNodesElement...
            (nodes,vertCoord,aNurbs,u1,u2);
%             nodesAdaptedRef = nefemInterp2DAdaptedNodesElement...
%                 (coordRef,vertCoord,aNurbs,u1,u2);
        nodesAdaptedRef = inverseLinearMapping(vertCoord,mesh.X(Te,:));

        % Vandermonde matrix
        V_nefem = Vandermonde_LP(nDeg,nodesAdaptedRef);
        invV_nefem = inv(V_nefem');

        % Compute shape functions at interpolation points
        shapeFunctions = zeros(npoints,nOfNodes);
        for ipoint = 1:npoints
            p = orthopoly2D(nodesAdapted(ipoint,:),nDeg);
            shapeFunctions(ipoint,:) = (invV_nefem*p)';
        end

        % Interpolate solution and position at interpolation points
        Xplot = shapeFunctions*mesh.X(Te,:);
        uplot = shapeFunctions*u(Te);

        % Plot interpolated solution in the element
        hold on
        elemPatchHandle(ielem) = patch('Faces',elemTriRef,'Vertices',Xplot,'FaceVertexCData',uplot,...
            'FaceColor','interp','EdgeAlpha',0);
        hold off
    end
    patchHandle = [patchHandle elemPatchHandle];
end
colorbarHandle = colorbar('location','East');
axis equal

%Output variable
if ~nargout
    varargout = [];
elseif nargout == 1
    varargout = {[patchHandle colorbarHandle]};
elseif nargout == 2
    varargout = {[patchHandle colorbarHandle] tri};
end
