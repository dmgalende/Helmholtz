function reflectingVector = berkhoffReflectingVector_CDG_P_adapt...
    (X,T,theReferenceElement,infoFaces,kvector,A,k,alpha,ABS,p_adaptivity)

elem_p = p_adaptivity.elem_p;
elem_check = p_adaptivity.elem_check;
allocation = 0.5*sum(((elem_p+1).*(elem_p+2)));
reflectingVector = zeros(allocation,1);
nOfElementNodes = size(theReferenceElement.NodesCoord,1);
nOfFaces = size(infoFaces,1);
face_check = elem_check(infoFaces(:,1));
faces = 1:nOfFaces;
% non refined element
for i = faces(~face_check)
    infoFace = infoFaces(i,:);
    iElem = infoFace(1);
    Xe = X(T(iElem,:),:);
    alpha_el = alpha(T(iElem,:));
    ke = k(T(iElem,:));
    ind = 0.5*sum((elem_p(1:iElem-1)+1).*(elem_p(1:iElem-1)+2))+1 : ...
        0.5*sum((elem_p(1:iElem)+1).*(elem_p(1:iElem)+2));
    rv_el = elementalReflectingVector(Xe,theReferenceElement,infoFace,kvector,A,ke,alpha_el,ABS,nOfElementNodes);
    reflectingVector(ind) = reflectingVector(ind) + rv_el;
end
% refined element
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
    rv_el = elementalReflectingVector(Xe,refEl,infoFace,kvector,A,ke,alpha_el,ABS,nOfElementNodes);
    reflectingVector(ind) = reflectingVector(ind) + rv_el;
end

function rv_el = elementalReflectingVector(Xe,theReferenceElement,infoFace,kvector,A,ke,ccge,ABS,nOfElementNodes)
iFace = infoFace(2);
nOfFaceNodes = numel(theReferenceElement.faceNodes1d);
rv = zeros(nOfFaceNodes,1);
rv_el = zeros(nOfElementNodes,1);
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
    %Incident potential at the current integration point
    xy_g = N_g*Xf;
    phi0_g = A*exp(sqrt(-1)*kvector*xy_g');
    phi0der_g = sqrt(-1)*kvector'*phi0_g;
    %Unit normal to the boundary
    t_g = xyDer_g/xyDerNorm_g;
    n_g = [t_g(2) -t_g(1)];
    %Contribution of the current integration point to the elemental vector
    rv = rv + ccg*N_g'*(n_g*phi0der_g - sqrt(-1)*k*ABS*phi0_g)*dline;
end
rv_el(nodes) = rv_el(nodes) + rv;