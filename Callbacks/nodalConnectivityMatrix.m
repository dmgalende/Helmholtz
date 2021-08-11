function N = nodalConnectivityMatrix(T)
 
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