function L2Norm = computeL2NormBoundary(theReferenceElement, X, T, u, uref_function, varargin)

% The function computeL2NormBoundaryScatteredCercle allows to compute the 
% L2 Norm on a boundary for a scattered wave in a cercle with a circular NRBC.
%
% Input:
%  referenceElement: information of the reference element
%  X: nodal coordinates
%  Tb: boundary connectivity matrix
%  u: FEM solution (nodal)ç
%  uref_function: handle to function computing the analytical solution f(rad, theta, opt_arg)
%  k: wave number
% Output:
%  L2Norm: computed L2 norm = sqrt(integral[(u)2]) over the boundary defined 
%          by Tb,X

%Number of elements and number of mesh nodes
nOfElements = size(T,1); 

%Loop in 2D elements
L2Norm = 0;
L2Norm_ref = 0;
for iElem = 1:nOfElements 
    Xe = X(T(iElem,:),:);
    ue = u(T(iElem,:),:);
    
    %Information of the reference element
    IPw = theReferenceElement.IPweights1d; 
    N = theReferenceElement.N1d; 
    Nxi = theReferenceElement.N1dxi;

    %Number of Gauss points
    ngauss = length(IPw);

    %Compute elemental L2 Norm
    for g = 1:ngauss
        %Values at current integration point
        N_g = N(g,:);
        Nxi_g = Nxi(g,:);
        ue_g = N_g*ue;
        %Analytical solution at integration point
        xy_g = N_g * Xe;
        [theta_g, rad_g] = cart2pol(xy_g(1), xy_g(2));
        uref_g = uref_function(rad_g, theta_g, varargin{:});
        %Integration weight
        xyDer_g = Nxi_g*Xe;
        xyDerNorm_g = norm(xyDer_g);
        dline=IPw(g)*xyDerNorm_g;
        %Contribution of the current integration point to the elemental L2 Norm
        err_g = ue_g - uref_g;
        L2Norm = L2Norm + err_g * conj(err_g) * dline;
        L2Norm_ref = L2Norm_ref + uref_g * conj(uref_g) * dline;
    end
    
end
L2Norm = sqrt(L2Norm / L2Norm_ref);
    
    
    
    
    
    
