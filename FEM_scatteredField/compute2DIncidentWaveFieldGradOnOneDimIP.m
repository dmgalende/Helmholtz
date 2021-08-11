function ugrad = compute2DIncidentWaveFieldGradOnOneDimIP(X,T,Tb,elementFaceInfo,...
                 referenceElement,kvector,amp0)

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
    
    %Choose shape functions
    if      iface == 1, N = N1;
    elseif  iface == 2, N = N2;
    else                N = N3;
    end
    
    %Loop in integration points
    for g = 1:ngauss
        
        %Shape function
        N_g = N(:,1,g)';

        %Grad of scalar field u at the current integration point
        xy_g = N_g*Xe;
        phi0_g = amp0*exp(sqrt(-1)*kvector*xy_g');
        phi0der_g = sqrt(-1)*kvector'*phi0_g;
        ugrad(g,1,ielem) = phi0der_g(1);
        ugrad(g,2,ielem) = phi0der_g(2);
        
    end
    
end


