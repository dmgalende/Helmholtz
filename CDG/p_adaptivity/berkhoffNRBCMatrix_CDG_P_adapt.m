function C = berkhoffNRBCMatrix_CDG_P_adapt(X,T,theReferenceElement,infoFaces,k,alpha,R,p_adaptivity)

nOfElementNodes = size(theReferenceElement.NodesCoord,1);
nOfFaces = size(infoFaces,1);
elem_p = p_adaptivity.elem_p;
elem_check = p_adaptivity.elem_check;
nOfElementNodes_vect = 0.5*((elem_p+1).*(elem_p+2));
allocation = sum(nOfElementNodes_vect(infoFaces(:,1)).^2);
dim = 0.5*sum((elem_p+1).*(elem_p+2));
I = zeros(allocation,1);
J = zeros(allocation,1);
Kr = zeros(allocation,1);
Ki = zeros(allocation,1);
aux_ones = ones(1,nOfElementNodes);
face_check = elem_check(infoFaces(:,1));
faces = 1:nOfFaces;
% non refined elements
for i = faces(~face_check)
    infoFace = infoFaces(i,:);
    iElem = infoFace(1);
    Xe = X(T(iElem,:),:);
    alpha_el = alpha(T(iElem,:));
    ke = k(T(iElem,:));
    ind = 0.5*sum((elem_p(1:iElem-1)+1).*(elem_p(1:iElem-1)+2))+1 : ...
        0.5*sum((elem_p(1:iElem)+1).*(elem_p(1:iElem)+2));
    Ce = elementalNRBCMatrix(Xe,theReferenceElement,infoFace,ke,alpha_el,R,nOfElementNodes);
    ind_transp = transpose(ind);
    aux_row = ind_transp(:,aux_ones); 
    aux_col = ind(aux_ones,:); 
    indK = 0.25*sum(((elem_p(1:iElem-1)+1).*(elem_p(1:iElem-1)+2)).^2)+1 : ...
        0.25*sum(((elem_p(1:iElem)+1).*(elem_p(1:iElem)+2)).^2);
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Kr(indK) = real(Ce); 
    Ki(indK) = imag(Ce); 
end
% refined elements
for i = faces(face_check)
    infoFace = infoFaces(i,:);
    iElem = infoFace(1);
    p = elem_p(iElem);
    refEl = p_adaptivity.referenceElements.(['P' num2str(p)]);
    nOfElementNodes = size(refEl.NodesCoord,1);
    indRef = 0.5*sum((elem_p(elem_check(1:iElem-1))+1).*(elem_p(elem_check(1:iElem-1))+2))+(1:nOfElementNodes);
    Xe = p_adaptivity.coord(indRef,:);
    alpha_el = p_adaptivity.ccg(indRef);
    ke = p_adaptivity.waveNumber(indRef);
    ind = 0.5*sum((elem_p(1:iElem-1)+1).*(elem_p(1:iElem-1)+2))+1 : ...
        0.5*sum((elem_p(1:iElem)+1).*(elem_p(1:iElem)+2));
    Ce = elementalNRBCMatrix(Xe,refEl,infoFace,ke,alpha_el,R,nOfElementNodes);
    ind_transp = transpose(ind);
    aux_row = ind_transp(:,ones(1,nOfElementNodes));
    aux_col = ind(ones(1,nOfElementNodes),:);
    indK = 0.25*sum(((elem_p(1:iElem-1)+1).*(elem_p(1:iElem-1)+2)).^2)+1 : ...
        0.25*sum(((elem_p(1:iElem)+1).*(elem_p(1:iElem)+2)).^2);
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Kr(indK) = real(Ce);
    Ki(indK) = imag(Ce);
end

C = sparse(I(I~=0),J(I~=0),Kr(I~=0)+sqrt(-1)*Ki(I~=0),dim,dim);


function Ce = elementalNRBCMatrix(Xe,theReferenceElement,infoFace,ke,ccge,R,nOfElementNodes)
iFace = infoFace(2);
nOfFaceNodes = numel(theReferenceElement.faceNodes1d);
C = zeros(nOfFaceNodes);
Ce = zeros(nOfElementNodes);
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
    %parameter for the Sommerfield NRB condition (radial boundary)
    param = sqrt(-1)*k - 1/(2*R);
    %Contribution of the current integration point to the elemental matrix
    C = C + ccg*param*N_g'*N_g*dline;
end
Ce(nodes,nodes) = Ce(nodes,nodes) + C; %elemental ensembling