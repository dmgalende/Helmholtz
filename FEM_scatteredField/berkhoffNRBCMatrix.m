function C = berkhoffNRBCMatrix(X,TRobin,theReferenceElement,k,ccg,sigma,omega,R)

% Input:
%  X: nodal coordinates
%  TRobin: connectivity matrix for 1D boundary elements
%  theReferenceElement: information of the reference element
%  k: nodal values of wave number
%  ccg: nodal values of the c*cg parameter
%  R: radius of circular non-reflecting boundary (radiation)
% Output:
%  C: non reflecting boundary condition matrix

%Number of elements and number of mesh nodes
nOfElements = size(TRobin,1);
nOfNodes = size(X,1);
nOfBoundaryNodes = length(unique(TRobin));
nOfLinearNodes = nOfElements + 1;
nOfElementNodes = theReferenceElement.degree + 1;
nOfInnerNodes = nOfBoundaryNodes - nOfLinearNodes;

%Memory allocation
NNZ = nOfElementNodes*(2*nOfLinearNodes + nOfInnerNodes);
% C = spalloc(nOfNodes,nOfNodes,NNZ);
aux_ones = ones(1,nOfElementNodes);
I = zeros(NNZ,1);
J = I;
Creal = I;
Cimag = I;
 
%Loop in 1D boundary elements
for ielem = 1:nOfElements
 Te = TRobin(ielem,:);
 Xe = X(Te,:); % Coordinates of the current element nodes
 ccge = ccg(Te);
 ke = k(Te);
 sigmae = sigma(Te,:);
 INDEX_PML_BOUNDARY = getIndex(Xe);
 Ce = ElementalNRBCMatrix(Xe,theReferenceElement,ke,ccge,sigmae,omega,R,INDEX_PML_BOUNDARY);
 
 % Assembling
 Te_transp = transpose(Te);
 aux_row = Te_transp(:,aux_ones);
 aux_col = Te(aux_ones,:);
 indK = (ielem-1)*nOfElementNodes^2+1:ielem*nOfElementNodes^2;
 I(indK) = aux_row(:);
 J(indK) = aux_col(:);
 Creal(indK) = real(Ce);
 Cimag(indK) = imag(Ce);
 
 clear Ce;
end

C = sparse(I(I~=0),J(I~=0),Creal(I~=0) + sqrt(-1)*Cimag(I~=0),nOfNodes,nOfNodes);

%______________________________________________________________
function Ce = ElementalNRBCMatrix(Xe,theReferenceElement,ke,ccge,sigmae,omega,R,INDEX_PML_BOUNDARY)

nOfElementNodes = size(Xe,1);
Ce = zeros(nOfElementNodes,nOfElementNodes);

%Information of the reference element
IPw = theReferenceElement.IPweights1d; 
N = theReferenceElement.N1d; 
Nxi = theReferenceElement.N1dxi;

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
  %phase and group celerities at current gauss point
  ccg = N_g*ccge;
  k = N_g*ke;
  %parameter for the Sommerfield NRB condition (radial boundary)
  param = sqrt(-1)*k - 1/(2*R);
  %PML parameter
  sigma_g = N_g*sigmae;
  param2 = 1 + sqrt(-1)*sigma_g(INDEX_PML_BOUNDARY)/omega;
  %Contribution of the current integration point to the elemental matrix
  Ce = Ce + ccg*param*param2*(N_g')*N_g*dline;
end


function i = getIndex(X)

tol = 1e-5;
j = 1;
res = tol/2;
while res < tol && j < size(X,1)
    res = abs((X(j,1) - X(j+1,1)) / X(j+1,1));
    j = j+1;
end
if res < tol, i = 2; else i = 1; end
   
    
