function L2Norm = computeL2NormBoundary(theReferenceElement,X,T,u)

% The function computeL2NormBoundaryScatteredCercle allows to compute the 
% L2 Norm on a boundary for a scattered wave in a cercle with a circular NRBC.
%
% Input:
%  referenceElement: information of the reference element
%  X: nodal coordinates
%  Tb: boundary connectivity matrix
%  u: FEM solution (nodal)
%  k: wave number
% Output:
%  L2Norm: computed L2 norm = sqrt(integral[(u)2]) over the boundary defined 
%          by Tb,X

%Number of elements and number of mesh nodes
nOfElements = size(T,1); 

%Loop in 2D elements
L2Norm = 0;
for iElem = 1:nOfElements 
    Xe = X(T(iElem,:),:);
    ue = u(T(iElem,:),:);
    L2Norm = L2Norm +  ElementalL2Norm(Xe,theReferenceElement,ue);  
end
L2Norm = sqrt(L2Norm);


%_______________________________________________________________________
function [elemL2Norm elemL2Norm_an] = ElementalL2Norm(Xe,theReferenceElement,ue)

%Information of the reference element
IPw = theReferenceElement.IPweights1d; 
N = theReferenceElement.N1d; 
Nxi = theReferenceElement.N1dxi;

%Number of Gauss points
ngauss = length(IPw);

%Compute elemental L2 Norm
elemL2Norm = 0;
elemL2Norm_an = 0;
for g = 1:ngauss
    %Values at current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    ue_g = N_g*ue;
    %Integration weight
    xyDer_g = Nxi_g*Xe;
    xyDerNorm_g = norm(xyDer_g);
    dline=IPw(g)*xyDerNorm_g;
    %Contribution of the current integration point to the elemental L2 Norm 
    elemL2Norm = elemL2Norm + abs(ue_g)^2*dline;
end
    
    
    
    
    
    
