function [A11,A12,A21,A22] = lifting_P_adapt(iFace_1,iFace_2,refEl_1,refEl_2,n,R_iElem,M_iFace,M_iFace_12,C12)

nOfElementNodes_1 = size(refEl_1.NodesCoord,1);
nOfElementNodes_2 = size(refEl_2.NodesCoord,1);
nodes = refEl_1.faceNodes(iFace_1,:);
%% calculate S11 
M_iFaceElem = zeros(nOfElementNodes_1);
M_iFaceElem(nodes,nodes) = M_iFace;
% % calculate M expanded
R_iElem_Expand = R_iElem;

% M_iFaceElem_Expand = expandMatrix(M_iFaceElem);
% the following lines replace expandMatrix 
nn = 2;
M_iFaceElem_Expand = zeros([size(M_iFaceElem) nn nn]);
M_iFaceElem_Expand(:,:,1:nn+1:nn^2) = repmat(M_iFaceElem, [1 1 nn]);
M_iFaceElem_Expand = permute(M_iFaceElem_Expand, [3 1 4 2]);
M_iFaceElem_Expand = reshape(M_iFaceElem_Expand, nn*size(M_iFaceElem));
%
S11 = zeros(2*nOfElementNodes_1,nOfElementNodes_1);
% apply the n1 normal
for j = 1: size(S11,2)
    col = j*2-1; % column index
    S11(:,j) = M_iFaceElem_Expand(:,col)*n(1) + M_iFaceElem_Expand(:,col+1)*n(2);
end

nodes2 = refEl_2.faceNodes(iFace_2,:);
%% calculate S12
M_iFace = fliplr(M_iFace_12);
M_iFaceElem = zeros(nOfElementNodes_1,nOfElementNodes_2);
M_iFaceElem(nodes,nodes2) =  M_iFace;
% M_iFaceElem_Expand = expandMatrix(M_iFaceElem);
% the following lines replace expandMatrix 
nn = 2;
M_iFaceElem_Expand = zeros([size(M_iFaceElem) nn nn]);
M_iFaceElem_Expand(:,:,1:nn+1:nn^2) = repmat(M_iFaceElem, [1 1 nn]);
M_iFaceElem_Expand = permute(M_iFaceElem_Expand, [3 1 4 2]);
M_iFaceElem_Expand = reshape(M_iFaceElem_Expand, nn*size(M_iFaceElem));
%
S12 = zeros(2*nOfElementNodes_1,nOfElementNodes_2);
% apply the n2 normal
for j = 1: size(S12,2)
    col = j*2-1; % column index
    S12(:,j) = -M_iFaceElem_Expand(:,col)*n(1) - M_iFaceElem_Expand(:,col+1)*n(2);
end
 %(n2 = -n1 )
A11 = transpose(S11)*R_iElem_Expand*S11 * (0.25+0.5*C12+(0.5*C12)^2);
A12 = transpose(S11)*R_iElem_Expand*S12 * (0.25+0.5*C12+(0.5*C12)^2);
A21 = transpose(A12);
A22 = transpose(S12)*R_iElem_Expand*S12 * (0.25+0.5*C12+(0.5*C12)^2);







