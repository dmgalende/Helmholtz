function varargout = plotSolutionInterpElem(mesh,u,computation,nDegRef,nameConRef)

% patchHandle = plotSolutionInterpElem(mesh,u,computation,nDegRef)
% plots a nodal scalar field using delaunay triangulation. The solution is
% interpolated on each element using a total number of interpolation points
% corresponding to an element of degree nDegRef (equal spaced points).
% 
% Input:
%   mesh: GUI mesh data structure.
%   u: nodal values.
%   computation: computation method (FEM, NEFEM, CDG, DG).
%   nDegRef (optional): degree of plotting element. By default is 15.
%   nameConRef (optional with nDegRef): name of the 2D mesh connectivity to be
%                                       plotted for FEM / CDG / DG computations
%
% Output:
%   patchHandle (optional): handle to the created patches objects

% Check input
if nargin == 3
    nDegRef = 15;
    nameCon = mesh.fieldNames{mesh.indexElemPosCon(3)};
elseif nargin == 4
    nameCon = mesh.fieldNames{mesh.indexElemPosCon(3)};
else
    nameCon = nameConRef;
end

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
%elemTriRef = delaunayn(nodes,{'QJ','QbB','Qc','Qx','Pp'});

% Vandermonde matrix
referenceElement = mesh.referenceElement;
coordRef = referenceElement.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = referenceElement.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');

% Select method and plot
if strcmp(computation,'FEM')
    try nameNodal = mesh.fieldNames{mesh.indexElemPosCon(2)};
    catch, nameNodal = 'X';
    end
    X = mesh.(nameNodal);
    T = mesh.(nameCon);
    [patchHandle, colorbarHandle] = plotInterpElemFEM(X,T);
elseif strcmp(computation,'NEFEM')
    [patchHandle, colorbarHandle] = plotInterpElemNEFEM();
elseif  any(strcmp(computation,{'CDG' 'DG'}))
    [patchHandle, colorbarHandle] = plotInterpElemCDG();
end

%Output variable
if ~nargout
    varargout = [];
elseif nargout == 1
    varargout = {[patchHandle colorbarHandle]};
end

%----------------------------

    function [patchHandle colorbarHandle] = plotInterpElemFEM(X,T)

    % Mesh info
    nOfElements = size(T,1);
    
    % Compute shape functions at interpolation points
    shapeFunctions = zeros(npoints,nOfNodes);
    for ipoint = 1:npoints
        p = orthopoly2D(nodes(ipoint,:),nDeg);
        shapeFunctions(ipoint,:) = (invV*p)';
    end
    
    % Loop in elements
    patchHandle = zeros(1,nOfElements);
    for ielem = 1:nOfElements

        % Interpolate solution and position at interpolation points
        Te = T(ielem,:);
        Xplot = shapeFunctions*X(Te,:);
        uplot = abs(shapeFunctions*real(u(Te)) + sqrt(-1)*shapeFunctions*imag(u(Te)));

        % Plot interpolated solution in the element
        hold on
        patchHandle(ielem) = patch('Faces',elemTriRef,'Vertices',Xplot,'FaceVertexCData',uplot,...
            'FaceColor','interp','EdgeAlpha',0);
        hold off
    end
    colorbarHandle = colorbar('location','East');
    axis equal

    end

%----------------------------

    function [patchHandle colorbarHandle] = plotInterpElemNEFEM()

    % Plot interior elements
    X = mesh.X;
    hold on
    [patchHandle colorbarHandle] = plotInterpElemFEM(X,mesh.T);
    hold off
    delete(colorbarHandle)
        
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
            vertCoord = X(iboundaryData(ielem,1:3),:);
            
            % Adapted reference element nodes to nurbs size
            nodesAdapted = nefemInterp2DAdaptedNodesElement...
                (nodes,vertCoord,aNurbs,u1,u2);
%             nodesAdaptedRef = nefemInterp2DAdaptedNodesElement...
%                 (coordRef,vertCoord,aNurbs,u1,u2);
            nodesAdaptedRef = inverseLinearMapping(vertCoord,X(Te,:));
            
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
            Xplot = shapeFunctions*X(Te,:);
            uplot = abs(shapeFunctions*real(u(Te)) + sqrt(-1)*shapeFunctions*imag(u(Te)));

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
    
    end

%----------------------------

    function [patchHandle colorbarHandle] = plotInterpElemCDG()

    % Mesh info
    nameNodal = mesh.fieldNames{mesh.indexElemPosCon(2)};
    nameCon = mesh.fieldNames{mesh.indexElemPosCon(3)};
    X = mesh.(nameNodal);
    T = mesh.(nameCon);
    nOfElementNodes = size(coordRef,1);
    nOfElements = size(T,1);
    
    % Compute shape functions at interpolation points
    shapeFunctions = zeros(npoints,nOfNodes);
    for ipoint = 1:npoints
        p = orthopoly2D(nodes(ipoint,:),nDeg);
        shapeFunctions(ipoint,:) = (invV*p)';
    end
    
    % Loop in elements
    patchHandle = zeros(1,nOfElements);
    for ielem = 1:nOfElements

        % Interpolate solution and position at interpolation points
        Te = T(ielem,:);
        Xplot = shapeFunctions*X(Te,:);
        ind = (ielem-1)*nOfElementNodes+1:ielem*nOfElementNodes;
        uplot = shapeFunctions*u(ind);

        % Plot interpolated solution in the element
        hold on
        patchHandle(ielem) = patch('Faces',elemTriRef,'Vertices',Xplot,'FaceVertexCData',uplot,...
            'FaceColor','interp','EdgeAlpha',0);
        hold off
    end
    colorbarHandle = colorbar('location','East');
    axis equal

    end

end
