function f = berkhoffExteriorVector(X,T,theReferenceElement,u0grad,ccg,sigma,omega)

% Input:
%  d0: unit vector of incident potential's direction
%  amp0: amplitude of incident potential
%  k: nodal values of wave number
%  ccg: nodal values of the c*cg parameter
%  X: nodal coordinates
%  T: connectivity matrix for 1D boundary elements
%  theReferenceElement: information of the reference element
%  alpha: parameter to set the reflection/absortion boundary condition
% Output:
%  f: berkhoff reflecting boundary vector

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfNodes = size(X,1);

%Memory allocation
f = zeros(nOfNodes,1);
 
%Loop in 1D boundary elements
for ielem = 1:nOfElements
 Te = T(ielem,:);
 Xe = X(Te,:); % Coordinates of the current element nodes
 ccge = ccg(Te);
 u0grade = u0grad(:,:,ielem);
 sigmae = sigma(Te,:);
 fe = forceElementalVector(Xe,theReferenceElement,u0grade,ccge,sigmae,omega);
 f(Te)=f(Te)+fe;
 clear fe;
end

%______________________________________________________________
function fe = forceElementalVector(Xe,theReferenceElement,u0grade,ccge,sigmae,omega)

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
  %Incident potential at the current integration point
  u0x_g = u0grade(g,1);
  u0y_g = u0grade(g,2);
  %PML parameters
  sigma_g = N_g*sigmae;
  param = 1 + sqrt(-1)*sigma_g/omega;
  %Unit normal to the boundary
  t_g = xyDer_g/xyDerNorm_g;
  nx_g = t_g(2);
  ny_g = -t_g(1);
  %Contribution of the current integration point to the elemental vector
  faux = ccg*((param(2)/param(1))*nx_g*u0x_g + (param(1)/param(2))*ny_g*u0y_g)*N_g';
  fe = fe + faux*dline;
end


