function C = coupledMassMatrix_u(X, Tb, T1d, refelem)

%Number of elements and number of mesh nodes
[nOfElements, nOfElementNodes] = size(Tb);
nOf1DElementNodes = size(T1d, 2);
nOf1DNodes = length(unique(T1d(:)));
nOfNodes = size(X, 1);

%Memory allocation
NNZ = nOfElementNodes * nOf1DElementNodes * nOfElements;
aux_ones_row = ones(1,nOf1DElementNodes);
aux_ones_col = ones(1,nOfElementNodes);
I = zeros(NNZ,1);
J = zeros(NNZ,1);
Creal = I;
Cimag = I;

%Loop in 1D boundary elements
for ielem = 1:nOfElements
    
    %Coords and numbering info
    Te = Tb(ielem,:);
    Te1d = T1d(ielem, :);
    Xe = X(Te, :);
    
    %Elemental matrix
    Ce = zeros(nOfElementNodes, nOf1DElementNodes);
    
    %Information of the reference element (assumed same basis for 2D and 1D)
    IPw = refelem.IPweights1d;
    N = refelem.N1d;
    Nxi = refelem.N1dxi;
    
    %Number of Gauss points
    ngauss = length(IPw);
    
    %Loop in integration points
    for g = 1:ngauss
        %Shape functions and derivatives at the current integration point
        N_g = N(g,:);
        Nxi_g = Nxi(g,:);
        %Integration weight
        xyDer_g = Nxi_g*Xe;
        xyDerNorm_g = norm(xyDer_g);
        dline=IPw(g)*xyDerNorm_g;
        %Contribution of the current integration point to the elemental matrix
        Ce = Ce + (N_g')*N_g*dline;
    end
    
    % Assembling
    Te_transp = transpose(Te);
    aux_row = Te_transp(:,aux_ones_row);
    aux_col = Te1d(aux_ones_col,:);
    indK = (ielem-1)*nOfElementNodes*nOf1DElementNodes+1:ielem*nOfElementNodes*nOf1DElementNodes;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Creal(indK) = real(Ce);
    Cimag(indK) = imag(Ce);
    
    clear Ce;
end

C = sparse(I(I~=0),J(I~=0),Creal(I~=0) + sqrt(-1)*Cimag(I~=0),nOfNodes,nOf1DNodes);
