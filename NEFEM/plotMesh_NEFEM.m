function varargout = plotMesh_NEFEM(mesh,npoints)

%
% patchHandle = plotMesh_NEFEM(mesh,npoints)
%
% Input
%   mesh: structure containing the X, T, Tb_, matrices and their 
%         corresponding indexes and name in the structure. Also the nurbs,  
%         the trimmedInfo and the reference element have to be inclouded.
%   npoints (optional): number of points the nurbs will be discretized with.
%                       By default are 20 points.
%
% Output
%   patchHandle (optional): handle (vector of handles) to the patch objects
%

% Dicretize the boundary face in npoints (always first face)
if nargin < 2
    npoints = 20;
end
faceNodes = mesh.referenceElement.faceNodes;
coordRef = mesh.referenceElement.NodesCoord;
xcoordRefDis = linspace(coordRef(faceNodes(1,1),1),coordRef(faceNodes(1,end),1),npoints);
ycoordRefDis = linspace(coordRef(faceNodes(1,1),2),coordRef(faceNodes(1,end),2),npoints);
coordRefDis = [xcoordRefDis' ycoordRefDis'];

% Nodes curved face (nurbs)
nodesNurbsRef = 1:npoints;

% Plot boundary elements
nurbs = mesh.nurbs;
nOfBoundaries = numel(mesh.boundaryNames);
patchHandle = zeros(1,nOfBoundaries+1);
for i = 1:nOfBoundaries
    item = mesh.boundaryIndex(i);
    iboundary = mesh.fieldNames{item};
    iboundaryTrimInfo = mesh.trimmedInfo.(iboundary);
    iboundaryData = mesh.(iboundary)(:,1:3);
    nOfElements = size(iboundaryData,1);
    
    % Allocate space for new nodes in nodal matrix
    extraSize = nOfElements*npoints;
    Xplot = [mesh.X ; zeros(extraSize,2)];
    ini = size(mesh.X,1) + 1;
    fin = ini + npoints - 1;
    
    % New connectivity (ordered by face)
    Tplot = zeros(nOfElements,npoints+1);
    Tplot(:,npoints+1) = iboundaryData(:,3); 
    
    for j = 1:nOfElements
        
        % NEFEM info
        u1 = iboundaryTrimInfo(j).trim(1);
        u2 = iboundaryTrimInfo(j).trim(2);
        idNurbs = iboundaryTrimInfo(j).idNurbs;
        aNurbs = nurbs(idNurbs);
        vertCoord = Xplot(iboundaryData(j,:),:);
        
        % Adapted reference element nodes (first face) to nurbs size
        nodesAdapted = nefemInterp2DAdaptedNodesElement...
            (coordRefDis,vertCoord,aNurbs,u1,u2);
        
        % Incloude nodes into new connectivity
        num = ini:fin;
        Tplot(j,nodesNurbsRef) = num;
        
        % Map to real space
        Xplot(num,:) = linearMapping(vertCoord,nodesAdapted);
        
        % Update nodes
        ini = fin + 1;
        fin = fin + npoints;
    end
    
    % Plot mesh (boundary elements)
    hold on
    patchHandle(i) = patch('Faces',Tplot,'Vertices',Xplot,'FaceColor',[1 1 1],'EdgeAlpha',1);
    hold off
end

% Plot mesh (interior elements with straight faces)
hold on
patchHandle(end) = patch('Faces',mesh.T(:,1:3),'Vertices',mesh.X,'FaceColor',[1 1 1],'EdgeAlpha',1);
hold off
axis equal

%Output variable
if ~nargout
    varargout = [];
else
    varargout = {patchHandle};
end