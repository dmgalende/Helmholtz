function C = coupledMassMatrix_f(X, T, T2d, refelem, periodic_elem, elementFaceInfo)

%Number of elements and number of mesh nodes
[nOfElements, nOfElementNodes] = size(T);
nOf2DElementNodes = size(T2d,2);
nOfNodes = length(unique(T(:)));
nOf2DNodes = length(unique(T2d(:)));

%Memory allocation
NNZ = nOfElementNodes * nOf2DElementNodes * nOfElements;
aux_ones_row = ones(1,nOf2DElementNodes);
aux_ones_col = ones(1,nOfElementNodes);
I = zeros(NNZ,1);
J = zeros(NNZ,1);
Creal = I;
Cimag = I;

%Loop in 1D boundary elements
for ielem = 1:nOfElements
    
    %Coords and numbering info
    Te = T(ielem,:);
    Te2d = T2d(elementFaceInfo(ielem,1), :);
    Xe = X(Te, :);
    
    %Fix angular coordinate in the element containing the periodic node
    if ielem == periodic_elem, Xe(end) = 2 * pi; end
    
    %Elemental matrix
    Ce = zeros(nOfElementNodes, nOf2DElementNodes);
    
    %Information of the reference element
    IPw = refelem.IPweights1d;
    N = refelem.N1d;
    Nxi = refelem.N1dxi;
    Nf = refelem.N_onFace1;
    
    %Number of Gauss points
    ngauss = length(IPw);
    
    %Loop in integration points
    for g = 1:ngauss
        %Shape functions and derivatives at the current integration point
        N_g = N(g,:);
        Nf_g = Nf(g,:);
        Nxi_g = Nxi(g,:);
        %Integration weight
        xyDer_g = Nxi_g*Xe;
        dline=IPw(g)*xyDer_g;
        %Contribution of the current integration point to the elemental matrix
        Ce = Ce + (N_g')*Nf_g*dline;
    end
    
    % Assembling
    Te_transp = transpose(Te);
    aux_row = Te_transp(:,aux_ones_row);
    aux_col = Te2d(aux_ones_col,:);
    indK = (ielem-1)*nOfElementNodes*nOf2DElementNodes+1:ielem*nOfElementNodes*nOf2DElementNodes;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Creal(indK) = real(Ce);
    Cimag(indK) = imag(Ce);
    
    clear Ce;
end

C = sparse(I(I~=0),J(I~=0),Creal(I~=0) + sqrt(-1)*Cimag(I~=0),nOfNodes,nOf2DNodes);
