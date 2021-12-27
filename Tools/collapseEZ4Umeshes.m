function collapsedMesh = collapseEZ4Umeshes(meshFile1,meshFile2,breakLineName,ortoVar)

% Input
if nargin == 3
   ortoVar = 1;
end

% Mesh parts
[~,~,ext1] = fileparts(meshFile1);
[~,~,ext2] = fileparts(meshFile2);
if all(strcmp({ext1 ext2},'.dcm'))
   matFiles = GenerateMatFileFromEZ4U({meshFile1 meshFile2});
   mesh1 = load(matFiles{1});
   mesh2 = load(matFiles{2});
elseif all(strcmp({ext1 ext2},'.mat'))
   mesh1 = load(meshFile1);
   mesh2 = load(meshFile2);
   matFiles = {[]};
elseif strcmp(ext1,'.mat')
   mesh1 = load(meshFile1);
   matFiles = GenerateMatFileFromEZ4U(meshFile2);
   mesh2 = load(matFiles{1});
else
   mesh2 = load(meshFile2);
   matFiles = GenerateMatFileFromEZ4U(meshFile1);
   mesh1 = load(matFiles{1});
end
delete(matFiles{:})

% Check doubled nodes
if max(max(mesh1.T)) ~= size(mesh1.X,1)
   error('The mesh 1 has doubled nodes!')
elseif max(max(mesh2.T)) ~= size(mesh2.X,1)
   error('The mesh 2 has doubled nodes!')
end

% Check the break line
tol = 1e-9;
breakLine = ['Tb_' breakLineName];
breakNodes1 = unique(mesh1.(breakLine));
breakNodes2 = unique(mesh2.(breakLine));
if length(breakNodes1) ~= length(breakNodes2)
   error('Number of break nodes in mesh 1 and mesh2 does not coincide')
end
[X1,pos1] = sort(mesh1.X(breakNodes1,ortoVar));
Y1 = mesh1.X(breakNodes1(pos1),3-ortoVar);
[X2,pos2] = sort(mesh2.X(breakNodes2,ortoVar));
Y2 = mesh2.X(breakNodes2(pos2),3-ortoVar);
xcondition = (X1 < X2 + tol) & (X1 > X2 - tol);
ycondition = (Y1 < Y2 + tol) & (Y1 > Y2 - tol);
breakCondition = xcondition & ycondition;
if ~all(breakCondition)
   disp(['The ' breakLine ' nodes...'])
   disp(pos1(~breakCondition))
   disp('do not coincide!!')
   error('Collapse not allowed')
else
   disp('Break line nodes coincide: collapse is allowed')
end

% Break nodes position in connectivity matrix
breakNodesValence = 5;
T2 = mesh2.T;
sizeT2 = size(T2);
nOfMesh2Nodes = size(mesh2.X,1);
N2 = getNodalConnectivity(T2);
breakNodes1 = breakNodes1(pos1);
breakNodes2 = breakNodes2(pos2);
nOfBreakNodes = length(breakNodes1);
breakNum = ones(nOfBreakNodes,1);
breakNodes2Pos = zeros(nOfBreakNodes,breakNodesValence);
for inode = 1:nOfBreakNodes
   ibreak2 = breakNodes2(inode);
   elem2 = N2(ibreak2,logical(N2(ibreak2,:)));
   for ielem = 1:length(elem2)
       breakElem = elem2(ielem);
       breakPos = find(T2(breakElem,:) == ibreak2);
       T2(breakElem,breakPos) = nOfMesh2Nodes + 1; %reset break nodes to maximum value
       breakNodes2Pos(inode,breakNum(inode)) = sub2ind(sizeT2,breakElem,breakPos);
       breakNum(inode) = breakNum(inode) + 1;
   end
end
breakNodes2Pos(:,max(breakNum):end) = [];

% New nodal position
X2 = mesh2.X;
X2(breakNodes2,:) = [];

% New connectivity
nOfMesh1Nodes = size(mesh1.X,1);
T2 = T2 + nOfMesh1Nodes; %with double break nodes
breakNodes2aux = sort(breakNodes2) + nOfMesh1Nodes;
for inode = 1:nOfBreakNodes
   ibreak2aux = breakNodes2aux(inode);
   nodesCondition = T2 > ibreak2aux;
   T2(nodesCondition) = T2(nodesCondition) - 1;
   breakNodes2aux(inode+1:nOfBreakNodes) = breakNodes2aux(inode+1:nOfBreakNodes) - 1;
end
minNum = min(min(T2));
T2 = T2 - minNum + nOfMesh1Nodes + 1; %without double break nodes

% Add break nodes of mesh1
for inode = 1:nOfBreakNodes
   ibreak1 = breakNodes1(inode);
   ibreak2Pos = breakNodes2Pos(inode,logical(breakNodes2Pos(inode,:)));
   T2(ibreak2Pos) = ibreak1;
end

% New mesh collapsed
mesh1 = rmfield(mesh1,['Tb_' breakLineName]);
mesh1.elementFaceInfo = rmfield(mesh1.elementFaceInfo,breakLineName);
mesh2 = rmfield(mesh2,['Tb_' breakLineName]);
mesh2.elementFaceInfo = rmfield(mesh2.elementFaceInfo,breakLineName);
collapsedMesh = mesh1;
collapsedMesh.T = [mesh1.T ; T2];
collapsedMesh.X = [mesh1.X ; X2];

% Boundaries
nOfMesh1Elements = size(mesh1.T,1);
faceNodes = mesh1.elemInfo.faceNodes;
fieldNames1 = fieldnames(mesh1);
nOfFieldNames1 = length(fieldNames1);
fieldNames2 = fieldnames(mesh2);
nOfFieldNames2 = length(fieldNames2);
for iboundary2 = 1:nOfFieldNames2
   iname2 = fieldNames2{iboundary2};
   if length(iname2) > 3 && strcmpi(iname2(1:3),'Tb_') && all(size(mesh2.(iname2)))
       Tb2 = mesh2.(iname2);
       % Renumbering the boundary connectivity
       elemFaceInfo2 = mesh2.elementFaceInfo.(iname2(4:end));
       nOfElements2 = size(Tb2,1);
       for ielem2 = 1:nOfElements2
           Tb2(ielem2,:) = T2(elemFaceInfo2(ielem2,1),faceNodes(elemFaceInfo2(ielem2,2),:));
       end
       elemFaceInfo2(:,1) = elemFaceInfo2(:,1) + nOfMesh1Elements;
       % Check for shared boundaries between mesh1 and mesh2
       iboundary1 = 1;
       iname1 = '';
       while ~strcmp(iname2,iname1) && iboundary1 <= nOfFieldNames1
           iname1 = fieldNames1{iboundary1};
           iboundary1 = iboundary1 + 1;
       end
       % Store in collapsedMesh       
       if strcmp(iname2,iname1) && all(size(mesh1.(iname1))) %shared
           collapsedMesh.(iname1) = [mesh1.(iname1) ; Tb2];
           collapsedMesh.elementFaceInfo.(iname1(4:end)) = ...
               [mesh1.elementFaceInfo.(iname1(4:end)) ; elemFaceInfo2];           
       else %not shared
           collapsedMesh.(iname2) = Tb2;
           collapsedMesh.elementFaceInfo.(iname2(4:end)) = elemFaceInfo2;
       end
   end
end


%---------------------------------------------------------
function N = getNodalConnectivity(T)

nNodes = max(max(T));
nNodesElem = size(T,2);
N = zeros(nNodes,10);
nn = ones(nNodes,1);
for ielem = 1:size(T,1)
   Te = T(ielem,:);
   nn_Te = nn(Te);
   for kk = 1:nNodesElem
       N(Te(kk),nn_Te(kk)) = ielem;
   end
   nn(Te) = nn(Te) + 1;
end
N(:,max(nn):end) = [];

           