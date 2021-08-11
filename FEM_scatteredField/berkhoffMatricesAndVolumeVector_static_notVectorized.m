function [K,f,elementalInfo,newPositionFaceNodes] = berkhoffMatricesAndVolumeVector_static...
    (X,T,theReferenceElement,amp0,kvector,k,ccg,sigma,omega,NNZ) 

% Input:
%  X: nodal coordinates
%  T: connectivity matrix for 2D elements
%  theReferenceElement: struct with the information of the reference
%  element
%  amp0: amplitude of incident potential
%  kvector: k0*d0, where k0, d0 are the wave number and the direction of
%  incident potential
%  k: nodal values of wave number
%  ccg: nodal values of the c*cg parameter
%  sigma: PML parameter
%  omega: angular frequency
%  NNZ: estimated number of non zero entries of the matrix
% Output:
%  K: (mass - stiffness) matrix ---ONLY FACE NODES---
%  f: volume vector ---ONLY FACE NODES---
%  elementalInfo: info needed for postprocess the inner nodes solution
%  newPositionFaceNodes: new global numbering for each face node which is 
%  going to be used to solve the system (zero numbering is used for inner
%  nodes)

%Reference element numbering
innerNodes = theReferenceElement.innerNodes;
nOfInnerNodes = length(innerNodes);
faceNodes = unique(theReferenceElement.faceNodes);
nOfFaceNodes = length(faceNodes);
  
%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfTotalInnerNodes = nOfElements*nOfInnerNodes;
nOfTotalNodes = size(X,1);
nOfNodes = nOfTotalNodes - nOfTotalInnerNodes;
%NNZ = NNZ - nOfTotalInnerNodes; %Fix estimated number of non zero entries
%if numel(NNZ) == 2 && NNZ(2) < 0, NNZ(2) = 0; end 

%Inner nodes position in global matrix (sorted)
innerNodesBool = false(1,nOfTotalNodes);
innerNodesBool(sort(reshape(T(:,innerNodes),nOfTotalInnerNodes,1),'ascend')) = true;

%New connectivity for face nodes (same sparsivity)
newPositionFaceNodes = cumsum(~innerNodesBool);
newPositionFaceNodes(innerNodesBool) = 0;

%Memory allocation
f = zeros(nOfNodes,1);
% K = spalloc(nOfNodes,nOfNodes,sum(NNZ));
% if numel(NNZ) == 2
%     Kreal = spalloc(nOfNodes,nOfNodes,NNZ(1));
%     Kimag = spalloc(nOfNodes,nOfNodes,NNZ(2));
% else
%     Kreal = spalloc(nOfNodes,nOfNodes,NNZ);
%     Kimag = spalloc(nOfNodes,nOfNodes,NNZ);
% end
nOfElementNodes = size(T,2);
nOfElementNodesMod = nOfElementNodes - nOfInnerNodes;
aux_ones = ones(1,nOfElementNodesMod);
allocation = nOfElementNodesMod^2*nOfElements;
I = zeros(allocation,1);
J = I;
Kreal = I;
Kimag = I;

%Elemental info
elementalInfo = zeros(nOfInnerNodes,nOfFaceNodes+1,nOfElements);

%Loop in 2D elements
for iElem = 1:nOfElements 
 Te = T(iElem,:); 
 Xe = X(Te,:);
 ccge = ccg(Te);
 ke = k(Te);
 sigmae = sigma(Te,:);
 [Ke,fe] = ElementalMatricesAndVolumeVectorBerkhoff...
     (Xe,theReferenceElement,amp0,kvector,ke,ccge,sigmae,omega);
 
 %Static condensation variables
 As = Ke(faceNodes,faceNodes);
 invBs = inv(Ke(innerNodes,innerNodes)); %inverse!
 cs = Ke(faceNodes,innerNodes);
 fs = fe(faceNodes);
 gs = fe(innerNodes);
 
 %New elemental matrix and vector (only face nodes)
 invBs_csT = invBs*cs.';
 invBs_gs = invBs*gs;
 Ke = As - cs*invBs_csT;
 fe = fs - cs*invBs_gs;
 
 %New connectivity for face nodes due to static condensation
 Te = newPositionFaceNodes(Te(faceNodes));
     
 %Store elemental info for postprocess the solution at inner nodes
 elementalInfo(:,:,iElem) = [invBs_gs invBs_csT];
 
 % Assembling
 f(Te) = f(Te) + fe;
%  K(Te,Te) = K(Te,Te) + Ke;
%  if any(any(sigmae))
%      SSO( Kimag, imag(Ke), Te, Te, '+' )
%      SSO( Kreal, real(Ke), Te, Te, '+' )
%  else
%      SSO( Kreal, Ke, Te, Te, '+' )
%  end
 Te_transp = transpose(Te);
 aux_row = Te_transp(:,aux_ones);
 aux_col = Te(aux_ones,:);
 indK = (iElem-1)*nOfElementNodesMod^2+1:iElem*nOfElementNodesMod^2;
 I(indK) = aux_row(:);
 J(indK) = aux_col(:);
 Kreal(indK) = real(Ke);
 Kimag(indK) = imag(Ke);
 
 clear Ke fe
end

% K = Kreal + sqrt(-1)*Kimag;
K = sparse(I,J,Kreal + sqrt(-1)*Kimag,nOfNodes,nOfNodes);


%__________________________________________________________________________
function [Ke,fe] = ElementalMatricesAndVolumeVectorBerkhoff...
    (Xe,theReferenceElement,amp0,kvector,ke,ccge,sigmae,omega) 
 
nOfElementNodes = size(Xe,1); %Number of nodes in the element

Ke = zeros(nOfElementNodes,nOfElementNodes); 
fe = zeros(nOfElementNodes,1);

%Information of the reference element
IPw = theReferenceElement.IPweights; 
N = theReferenceElement.N; 
Nxi = theReferenceElement.Nxi;
Neta = theReferenceElement.Neta;

%Number of Gauss points
ngauss = length(IPw); 
 
% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

%LOOP IN GAUSS POINTS
for g = 1:ngauss
  %Shape functions and derivatives at the current integration point 
  N_g = N(g,:);  
  Nxi_g = Nxi(g,:);   
  Neta_g = Neta(g,:);
  %Jacobian
  J = [Nxi_g*xe	  Nxi_g*ye   
       Neta_g*xe  Neta_g*ye];
  %Integration weight
  detJ = det(J);
  dvolu=IPw(g)*detJ;
  %x and y derivatives
  invJ11 = J(2,2); %Not 1/detJ (optimized)!
  invJ12 = -J(1,2);
  invJ21 = -J(2,1);
  invJ22 = J(1,1);
  invdetJ = 1/detJ;
  Nx_g = invdetJ*(invJ11*Nxi_g + invJ12*Neta_g);
  Ny_g = invdetJ*(invJ21*Nxi_g + invJ22*Neta_g);
  %phase and group celerities at current gauss point
  ccg = N_g*ccge;
  ccgx = Nx_g*ccge;
  ccgy = Ny_g*ccge;
  ccgGrad = [ccgx ccgy];
  k = N_g*ke;
  %Incident potential and derivative at the current integration point
  xy_g = N_g*Xe;
  phi0_g = amp0*exp(1i*kvector*xy_g');
  phi0der_g = sqrt(-1)*kvector'*phi0_g;
  phi0Lap_g = -kvector*kvector'*phi0_g;
  %PML parameters
  sigma_g = N_g*sigmae;
  param = 1 + sqrt(-1)*sigma_g/omega;
  %Contribution of the current integration point to the elemental matrix
  K_e = (dvolu*ccg)*((param(2)/param(1))*Nx_g'*Nx_g + (param(1)/param(2))*Ny_g'*Ny_g);
  M_e = (ccg*k^2*param(1)*param(2)*dvolu)*(N_g'*N_g);
  Ke = Ke + M_e - K_e;
  fe = fe + (N_g'*((param(1)*param(2)*dvolu)*(ccgGrad*phi0der_g + ccg*phi0Lap_g + ...
      k^2*ccg*phi0_g)));
end

 
