function [OutOfPMLelements,OutOfPMLelements_contact,N] = createOutOfPMLConnectivity(data)

X = data.mesh.X;
T = data.mesh.T;
nOfBoundaries = numel(data.mesh.boundaryNames);
numbering = 1:size(X,1);
elements = 1:size(T,1);

% n of PMLs
totalPML = [];
for iboundary = 1:nOfBoundaries
   if strcmp(data.PML{5,iboundary},'on')
       totalPML = [totalPML iboundary];
   end
end
nOfPMLs = numel(totalPML);

% nodal connectivity
nNodes = max(max(T));
N = zeros(nNodes,10);
nn = ones(nNodes,1);
for ielem = 1:size(T,1)
   Te = T(ielem,:);
   nn_Te = nn(Te);
   for kk = 1:3
       N(Te(kk),nn_Te(kk)) = ielem;
   end
   nn(Te) = nn(Te) +1;
end
N(:,max(nn):end) = [];

%PML nodes
nodesPML = false(size(X,1),1);
for i = 1:nOfPMLs
   nodesPMLaux = false(size(X,1),1);
   param = data.BC.parameters{totalPML(i)}{4}(4);
   if param > 1
       nameSubCon = data.mesh.fieldNames{data.mesh.submeshIndex(param-1)};
       nodesPMLaux(getNodes(data.mesh.(nameSubCon))) = true;
   else
       nodesPMLaux(data.PML{4,totalPML(i)}{1}) = true;
       nodesPMLaux(data.PML{4,totalPML(i)}{2}) = true;
   end
   nodesPML = nodesPML | nodesPMLaux;
end

% Out of PML elements
nodesPMLnum = numbering(nodesPML);
nOfElements = size(T,1);
nOfTimesMarked = zeros(nOfElements,4);
for inode = nodesPMLnum
   elem = N(inode,logical(N(inode,:)));
   for j = 1:length(elem)
       jelem = elem(j);
       jpos = find(T(jelem,:) == inode);
       nOfTimesMarked(jelem,1) = nOfTimesMarked(jelem,1) + 1;
       nOfTimesMarked(jelem,jpos+1) = 1;
   end
end
inPML = nOfTimesMarked(:,1) == 3; %Only elements with all PML nodes
onPML = nOfTimesMarked(:,1) == 2; %Only elements with 2 PML nodes
OutOfPMLelements = elements(~inPML);
OutOfPMLelements_contact = elements(onPML);
