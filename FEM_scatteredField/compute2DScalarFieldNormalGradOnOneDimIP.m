function ugradN = compute2DScalarFieldNormalGradOnOneDimIP(u,X,T,Tb,elementFaceInfo,referenceElement,order)

%Set the normal direction at the boundary defined by Tb
switch order
    case 'normal'
        nfactor = -1;
    case 'reverse'
        nfactor = 1;
end

%General
nOfElements = size(Tb,1);
nOfNodesElement = size(T,2);

%Information of the reference element
IPw = referenceElement.IPweights1d;
nDeg = referenceElement.degree;
ngauss = length(IPw);
Nxi1D = theReferenceElement.N1dxi;

%Integration points for 2d element (depend on the face)
ipcoordinates1 = one2twoDmapping(1,referenceElement);
ipcoordinates2 = one2twoDmapping(2,referenceElement);
ipcoordinates3 = one2twoDmapping(3,referenceElement);

%2D shape functions on the 1D integration points (depend on the face)
V = Vandermonde_LP(nDeg,referenceElement.NodesCoord);
[L,U,P] = lu(V');
N1 = zeros(nOfNodesElement,2,ngauss);
N2 = N1;
N3 = N1;
for g = 1:ngauss
    [~,p_xi,p_eta] = orthopoly2D_deriv_xieta(ipcoordinates1(g,:),nDeg);
    N1(:,:,g) = U\(L\(P*[p_xi,p_eta]));
    [~,p_xi,p_eta] = orthopoly2D_deriv_xieta(ipcoordinates2(g,:),nDeg);
    N2(:,:,g) = U\(L\(P*[p_xi,p_eta]));
    [~,p_xi,p_eta] = orthopoly2D_deriv_xieta(ipcoordinates3(g,:),nDeg);
    N3(:,:,g) = U\(L\(P*[p_xi,p_eta]));
end

%Loop in 1D boundary elements
ugradN = zeros(ngauss,nOfElements);
for ielem = 1:nOfElements
    
    %Element info
    iface = elementFaceInfo(ielem,2);
    jelem = elementFaceInfo(ielem,1);
    Te = T(jelem,:); %Nodes of the 2D element nodes
    Xe = X(Te,:);    %Coordinates of the 2D element nodes
    ue = u(Te);
    Xe1D = X(Tb(ielem,:),:); %Coordinates of the 1D element nodes
    
    %Choose shape functions
    if      iface == 1, N = N1;
    elseif  iface == 2, N = N2;
    else                N = N3;
    end
    
    %x and y coordinates of the 2D element
    xe = Xe(:,1); ye = Xe(:,2);
    
    %Loop in integration points
    for g = 1:ngauss
        
        %1D shape function
        Nxi_g = Nxi1D(g,:);
        
        %Jacobian
        N2xi_g = N(:,1,g)';
        N2eta_g = N(:,2,g)';
        J = [N2xi_g*xe	  N2xi_g*ye
            N2eta_g*xe  N2eta_g*ye];
        detJ = det(J);
        
        %x and y derivatives
        invJ11 = J(2,2);
        invJ12 = -J(1,2);
        invJ21 = -J(2,1);
        invJ22 = J(1,1);
        invdetJ = 1/detJ;
        Nx_g = invdetJ*(invJ11*N2xi_g + invJ12*N2eta_g);
        Ny_g = invdetJ*(invJ21*N2xi_g + invJ22*N2eta_g);
        
        %Unit normal to the boundary
        xyDer_g = Nxi_g*Xe1D;
        xyDerNorm_g = norm(xyDer_g);
        t_g = xyDer_g/xyDerNorm_g;
        n_g = [t_g(2) nfactor*t_g(1)];

        %Normal variation of scalar field u at the current integration point
        ugrad = [Nx_g*ue ; Ny_g*ue];
        ugradN(g,ielem) = n_g*ugrad; 
    end
    
end


