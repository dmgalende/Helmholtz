function [K,M] = matrices1D_Berkhoff(X,T,theReferenceElement,sigma,ccg,k,kx,omega)

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
M = I;
K = I;

%Information of the reference element
IPw = theReferenceElement.IPweights1d; 
N = theReferenceElement.N1d; 
Nxi = theReferenceElement.N1dxi;

%Number of Gauss points
ngauss = length(IPw);

%Loop in 1D boundary elements
for ielem = 1:nOfElements
    Te = T(ielem,:);
    Xe = X(Te,:);
    ke = k(Te);
    ccge = ccg(Te);
    sigmae = sigma(Te);
    
    %Loop in gauss points
    Me = zeros(nOfElementNodes,nOfElementNodes);
    Ke = Me;
    for g = 1:ngauss
        %Shape functions and derivatives at the current integration point
        N_g = N(g,:);
        Nxi_g = Nxi(g,:);
        k_g = N_g*ke;
        ccg_g = N_g*ccge;
        %Integration weight
        xyDer_g = Nxi_g*Xe;
        xyDerNorm_g = norm(xyDer_g);
        dline=IPw(g)*xyDerNorm_g;
        %PML parameters
        sigma_g = N_g*sigmae;
        param_u = 1 + sqrt(-1)*sigma_g/omega;
        %Contribution of the current integration point to the elemental matrix
        Maux = (N_g')*N_g*dline;
        Me = Me + param_u*ccg_g*(k_g^2-kx^2)*Maux;
        Kaux = (1/(xyDer_g^2))*(Nxi_g')*Nxi_g*dline;
        Ke = Ke + (1/param_u)*ccg_g*Kaux;
    end
    
    % Assembling
    Te_transp = transpose(Te);
    aux_row = Te_transp(:,aux_ones);
    aux_col = Te(aux_ones,:);
    indK = (ielem-1)*nOfElementNodes^2+1:ielem*nOfElementNodes^2;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    M(indK) = Me;
    K(indK) = Ke;
    
    clear Me Ke
end

%Create sparse output matrices
pos = I~=0;
M = sparse(I(pos),J(pos),M(pos),nOfNodes,nOfNodes);
K = sparse(I(pos),J(pos),K(pos),nOfNodes,nOfNodes);



