function f = berkhoffVolumeVector_notVectorized...
    (X,T,theReferenceElement,u0,kVec,ccgVec)

% Input:
%  X: nodal coordinates
%  T: connectivity matrix for 2D elements
%  theReferenceElement: struct with the information of the reference
%  element
%  amp0: amplitude of incident potential
%  kvector: k0*d0, where k0, d0 are the wave number and the direction of
%  incident potential
%  kVec: nodal values of wave number
%  ccgVec: nodal values of the c*cg parameter
%  sigmaVec: PML parameter
%  omega: angular frequency
%  NNZ: estimated number of non zero entries of the matrix
% Output:
%  K: (mass - stiffness) matrix
%  f: volume vector

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfNodes = size(X,1);

%Memory allocation
f = zeros(nOfNodes,1);

%Loop in 2D elements
for iElem = 1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    ccge = ccgVec(Te);
    ke = kVec(Te);
    u0e = u0(Te);
    
    fe = ElementalVolumeVectorBerkhoff(Xe,theReferenceElement,u0e,ke,ccge);

    % Assembling
    f(Te) = f(Te) + fe;
    
    clear fe
end


%__________________________________________________________________________
function fe = ElementalVolumeVectorBerkhoff(Xe,theReferenceElement,u0e,ke,ccge) 
 
nOfElementNodes = size(Xe,1); %Number of nodes in the element

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
  k = N_g*ke;
  %Incident potential and derivative at the current integration point
  u0_g = N_g*u0e;
  u0x_g = Nx_g*u0e;
  u0y_g = Ny_g*u0e;
  %Contribution of the current integration point to the elemental matrix
  faux = ccg*((k^2*u0_g)*N_g' - (u0x_g*Nx_g' + u0y_g*Ny_g'));
  fe = fe + faux*dvolu;
end
