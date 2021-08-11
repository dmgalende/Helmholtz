function K = stiffMatrix1D(X, T, theReferenceElement, periodic_elem)

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfNodes = size(X,1);
nOfBoundaryNodes = length(unique(T));
nOfLinearNodes = nOfElements + 1;
nOfElementNodes = theReferenceElement.degree + 1;
nOfInnerNodes = nOfBoundaryNodes - nOfLinearNodes;

%Memory allocation
NNZ = nOfElementNodes*(2*nOfLinearNodes + nOfInnerNodes);
aux_ones = ones(1,nOfElementNodes);
I = zeros(NNZ,1);
J = I;
K = I;

%Information of the reference element
IPw = theReferenceElement.IPweights1d;
Nxi = theReferenceElement.N1dxi;

%Number of Gauss points
ngauss = length(IPw);

%Loop in 1D boundary elements
for ielem = 1:nOfElements
    Te = T(ielem,:);
    Xe = X(Te,:);
    
    %Fix angular coordinate in the element containing the periodic node
    if ielem == periodic_elem, Xe(end) = 2 * pi; end
    
    %Loop in gauss points
    Ke = zeros(nOfElementNodes,nOfElementNodes);
    for g = 1:ngauss
        %Shape functions and derivatives at the current integration point
        Nxi_g = Nxi(g,:);
        %Integration weight
        xyDer_g = Nxi_g*Xe;
        dline=IPw(g)*xyDer_g;
        %Contribution of the current integration point to the elemental matrix
        Kaux = (1/xyDer_g^2) * (Nxi_g')*Nxi_g*dline;
        Ke = Ke + Kaux;
    end
    
    % Assembling
    Te_transp = transpose(Te);
    aux_row = Te_transp(:,aux_ones);
    aux_col = Te(aux_ones,:);
    indK = (ielem-1)*nOfElementNodes^2+1:ielem*nOfElementNodes^2;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    K(indK) = Ke;
    
    clear Ke
end

%Create sparse output matrices
pos = I~=0;
K = sparse(I(pos),J(pos),K(pos),nOfNodes,nOfNodes);



