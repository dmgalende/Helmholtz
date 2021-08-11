function dampingMatrix = berkhoffDampingMatrix_CDG(X,T,theReferenceElement,infoFaces,k,alpha)

nOfElementNodes = size(theReferenceElement.NodesCoord,1);
nOfFaces = size(infoFaces,1);
allocation = nOfFaces*nOfElementNodes;
I = zeros(allocation,1);
J = zeros(allocation,1);
Kr = zeros(allocation,1); 
Ki = zeros(allocation,1);
aux_ones = ones(1,nOfElementNodes);
for i = 1:nOfFaces
    infoFace = infoFaces(i,:);
    iElem = infoFace(1);
    Xe = X(T(iElem,:),:);
    alpha_el = alpha(T(iElem,:));
    ke = k(T(iElem,:));
    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    dm_el = elementalDampinMatrix(Xe,theReferenceElement,infoFace,ke,alpha_el,nOfElementNodes);
    ind_transp = transpose(ind);
    aux_row = ind_transp(:,aux_ones); 
    aux_col = ind(aux_ones,:); 
    indK = (iElem-1)*nOfElementNodes^2+1:iElem*nOfElementNodes^2;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Kr(indK) = real(dm_el); 
    Ki(indK) = imag(dm_el); 
end
dampingMatrix = sparse(I(I~=0),J(I~=0),Kr(I~=0)+sqrt(-1)*Ki(I~=0),numel(T),numel(T));


function dm_el = elementalDampinMatrix(Xe,theReferenceElement,infoFace,ke,ccge,nOfElementNodes)
iFace = infoFace(2);
nOfFaceNodes = numel(theReferenceElement.faceNodes1d);
dm = zeros(nOfFaceNodes);
dm_el = zeros(nOfElementNodes);
%Information of the reference element
IPw = theReferenceElement.IPweights1d; 
N = theReferenceElement.N1d; 
Nxi = theReferenceElement.N1dxi;

nodes = theReferenceElement.faceNodes(iFace,:);
Xf = Xe(nodes,:);
kf = ke(nodes,:);
ccgf = ccge(nodes,:);
%Number of Gauss points
ngauss = length(IPw);
%Loop in integration points
for g = 1:ngauss
    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    %Integration weight
    xyDer_g = Nxi_g*Xf;
    xyDerNorm_g = norm(xyDer_g);
    dline=IPw(g)*xyDerNorm_g;
    %phase and group celerities at current gauss point
    ccg = N_g*ccgf;
    k = N_g*kf;
    %Contribution of the current integration point to the elemental matrix
    dm = dm + ccg*k*N_g'*N_g*dline;
end
dm_el(nodes,nodes) = dm_el(nodes,nodes) + dm;