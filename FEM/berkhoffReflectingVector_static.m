function f = berkhoffReflectingVector_static(X,T,theReferenceElement,kvector,amp0,k,ccg,alpha,faceNodesNum,nOfNodes)

% Input:
%  d0: unit vector of incident potential's direction
%  amp0: amplitude of incident potential
%  k: nodal values of wave number
%  ccg: nodal values of the c*cg parameter
%  X: nodal coordinates
%  T: connectivity matrix for 1D boundary elements
%  theReferenceElement: information of the reference element
%  alpha: parameter to set the reflection/absortion boundary condition
%  faceNodesNum: new face nodes numbering used to solve the system
%  nOfNodes: allocation space accounting for static condensation
% Output:
%  f: berkhoff reflecting boundary vector

%Number of elements
nOfElements = size(T,1);

%Memory allocation
f = zeros(nOfNodes,1);
 
%Loop in 1D boundary elements
for ielem = 1:nOfElements
 Te = T(ielem,:);
 Xe = X(Te,:); % Coordinates of the current element nodes
 ccge = ccg(Te);
 ke = k(Te);
 fe = forceElementalVector(Xe,theReferenceElement,kvector,amp0,ke,ccge,alpha);
 
 %New face nodes numbering due to static condensation
 Te = faceNodesNum(Te);
 
 f(Te)=f(Te)+fe;
 clear fe;
end

%______________________________________________________________
function fe = forceElementalVector(Xe,theReferenceElement,kvector,amp0,ke,ccge,alpha)

nOfElementNodes = size(Xe,1);
fe = zeros(nOfElementNodes,1);

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
  %Incident potential at the current integration point
  xy_g = N_g*Xe;
  phi0_g = amp0*exp(sqrt(-1)*kvector*xy_g');
  phi0der_g = sqrt(-1)*kvector'*phi0_g;
  %Unit normal to the boundary
  t_g = xyDer_g/xyDerNorm_g;
  n_g = [t_g(2) -t_g(1)];
  %Contribution of the current integration point to the elemental vector
  fe = fe + ccg*N_g'*(n_g*phi0der_g - sqrt(-1)*k*alpha*phi0_g)*dline;
end


