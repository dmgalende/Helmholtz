function v = getNodes(T)

% Input:
%  T: connectivity matrix of the mesh
% Output:
%  v: vector containing the ordered mesh nodes without connectivity

%Get matrix information
[nOfElements,nOfNodesPerElement] = size(T);
nOfTotalNodes = nOfElements*nOfNodesPerElement;
v_ = zeros(1,nOfTotalNodes);

%Create 1-row vector
aux = [1 nOfNodesPerElement];
for ielem = 1:nOfElements
    inodes = T(ielem,:);
    v_(aux(1):aux(2)) = inodes;
    aux = aux + nOfNodesPerElement;  
end

%Order the 1-row vector and eliminate connectivities
v = unique(v_);
    
    
        