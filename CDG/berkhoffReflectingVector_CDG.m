function reflectingVector = berkhoffReflectingVector_CDG(X,T,theReferenceElement,infoFaces,kvector,A,k,alpha,ABS)

reflectingVector = zeros(numel(T),1);
nOfElementNodes = size(theReferenceElement.NodesCoord,1);
nOfExteriorPEC = size(infoFaces,1);
for i = 1:nOfExteriorPEC
    infoFace = infoFaces(i,:);
    iElem = infoFace(1);
    Xe = X(T(iElem,:),:);
    alpha_el = alpha(T(iElem,:));
    ke = k(T(iElem,:));
    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    rv_el = elementalReflectingVector(Xe,theReferenceElement,infoFace,kvector,A,ke,alpha_el,ABS,nOfElementNodes);
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