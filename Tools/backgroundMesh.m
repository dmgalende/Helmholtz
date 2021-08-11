
%% LOAD meshes

% boundaryMeshPath = '../../Harbors_PML/elliptic_shoal/geo/';
boundaryMeshPath = '../';
boundaryMeshFile = 'BoundaryMesh_NEFEM1.dcm';
bottomData       = '../BoundaryMesh_NEFEM1.mat';
matfiles = GenerateMatFileFromEZ4U(boundaryMeshPath,boundaryMeshFile);
load(matfiles{1})
delete(matfiles{:})

%% INTERPOLATION of bottom

%Initialization
d = nan(size(X,1),1);

%Interpolate bottom
load(bottomData)
D = scatteredInterpolant(data.mesh.X(:,1),data.mesh.X(:,2),data.bottom.value,'linear','nearest');
domainNodes = unique(T);
d(domainNodes) = D(X(domainNodes,1),X(domainNodes,2));
if any(d<0), warning('Negative values of bottom interpolation detected'), d(d<0) = min(d(d>0)); end

%Correction in the PML if needed
% pmlNodes = unique(data.mesh.extT);
% d(pmlNodes) = max(d(domainNodes));

%Boundary
Tboundary = [];
boundaryNames = fieldnames(elementFaceInfo);
for i = 1:length(boundaryNames)
    Tboundary = [Tboundary ; eval(['Tb_' boundaryNames{i}])];
end

%% COMPUTE ELEMENT SIZE

%Element size given by the bottom
period = 0.63*60;
N = 10;
p = 4;
L = 2*pi./computeWaveNumber(2*pi/period,d);
h = L * p / N;

%Element size fixed on the boundary
tol = 1e-3;
correction = 0.7;
nElemsBoundary = size(Tboundary,1);
N = nodalConnectivityMatrix(T);
nodesBoundary = unique(Tboundary);
nNodesBoundary = numel(nodesBoundary);
for i = 1:nNodesBoundary
    inode = nodesBoundary(i);
    ielems = N(inode,logical(N(inode,:)));
    Te = T(ielems,:);
    neighbours = setdiff(Te(:),inode);
    Xneighbours = X(neighbours,:);
    Xnode = X(inode,:);
    Xdiff = Xnode(ones(numel(neighbours),1),:) - Xneighbours;
    dis = sqrt(sum(Xdiff.^2,2)) * correction;
    h(inode) = min([h(inode) ; dis]);
end

%Correction outside the domain with the maximum value
% h(isnan(h)) = max(h(~isnan(h)));

%Delaunay triangulation without constraints for the best quality background mesh
DT = delaunayTriangulation(X(:,1),X(:,2));
DTi = scatteredInterpolant(X(:,1),X(:,2),h,'linear','nearest');
hi = DTi(DT.Points(:,1),DT.Points(:,2));

%%  WRITE BACKGROUND MESH

nNodes = size(DT.Points,1);
nElems = size(DT.ConnectivityList,1);
nodeId = 1:nNodes;
elemId = 1:nElems;
fid = fopen('backgroundMesh.sub','w');
fprintf(fid,'%i\n%i\n',nNodes,nElems);
fprintf(fid,'%i %g %g %g\n',[nodeId ; DT.Points' ; hi']);
fprintf(fid,'%i %i %i %i\n',[elemId ; DT.ConnectivityList']);
fclose(fid);






