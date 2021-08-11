function ugrad = compute2DScalarFieldGradOnOneDimIP(u,X,T,Tb,elementFaceInfo,referenceElement)

%General
nOfElements = size(Tb,1);
nOfNodesElement = size(T,2);

%Information of the reference element
IPw = referenceElement.IPweights1d;
nDeg = referenceElement.degree;
ngauss = length(IPw);

%Integration points for 2d element (depend on the face)
ipcoordinates1 = one2twoDmapping(1,referenceElement);
ipcoordinates2 = one2twoDmapping(2,referenceElement);
ipcoordinates3 = one2twoDmapping(3,referenceElement);

%2D shape functions on the 1D integration points (depend on the face)
V = Vandermonde_LP(nDeg,referenceElement.NodesCoord);
[L,U,P] = lu(V');
N1 = zeros(nOfNodesElement,3,ngauss);
N2 = N1;
N3 = N1;
for g = 1:ngauss
    [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(ipcoordinates1(g,:),nDeg);
    N1(:,:,g) = U\(L\(P*[p,p_xi,p_eta]));
    [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(ipcoordinates2(g,:),nDeg);
    N2(:,:,g) = U\(L\(P*[p,p_xi,p_eta]));
    [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(ipcoordinates3(g,:),nDeg);
    N3(:,:,g) = U\(L\(P*[p,p_xi,p_eta]));
end

%Loop in 1D boundary elements
ugrad = zeros(ngauss,2,nOfElements);
for ielem = 1:nOfElements
    
    %Element info
    iface = elementFaceInfo(ielem,2);
    jelem = elementFaceInfo(ielem,1);
    Te = T(jelem,:); %Nodes of the 2D element nodes
    Xe = X(Te,:);    %Coordinates of the 2D element nodes
    ue = u(Te);
    
    %Choose shape functions
    if      iface == 1, N = N1;
    elseif  iface == 2, N = N2;
    else                N = N3;
    end
    
    %x and y coordinates of the 2D element
    xe = Xe(:,1); ye = Xe(:,2);
    
    %Loop in integration points
    for g = 1:ngauss
        
        %Jacobian
        N2xi_g = N(:,2,g)';
        N2eta_g = N(:,3,g)';       
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

        %Grad of scalar field u at the current integration point
        ugrad(g,1,ielem) = Nx_g*ue;
        ugrad(g,2,ielem) = Ny_g*ue;
    end
    
end


