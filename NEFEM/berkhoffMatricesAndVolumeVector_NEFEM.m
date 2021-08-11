function [K,f] = berkhoffMatricesAndVolumeVector_NEFEM...
    (X,T,referenceElement,nurbs,trimmedInfo,amp0,kvector,k,ccg,sigma,omega,NNZ_nefem)

%
% [K,f] = berkhoffMatricesAndVolumeVector_NEFEM...
%    (X,T,referenceElement,nurbs,trimmedInfo,amp0,k,ccg,sigma,omega)
%
% Output K: (mass - stiffness) matrix
%

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfNodes = size(X,1);
nDeg = referenceElement.degree;

%Memory allocation
f = zeros(nOfNodes,1);
nOfElementNodes = size(T,2);
aux_ones = ones(1,nOfElementNodes);
allocation = nOfElementNodes^2*nOfElements;
I = zeros(allocation,1);
J = I;
Kreal = I;
Kimag = I;

% Gauss quadrature options
global nIPNurbsEdge nIPinterior

% FEM reference element
coordRef = referenceElement.NodesCoord;

% Vandermonde matrix
%     V = Vandermonde_LP(nDeg,coordRef);
%     invV = inv(V');

%Loop in 2D elements
for iElem = 1:nOfElements
    Te = T(iElem,:);
    Xe = X(Te,:);
    ccge = ccg(Te);
    ke = k(Te);
    sigmae = sigma(Te,:);

    % NEFEM information
    u1 = trimmedInfo(iElem).trim(1);
    u2 = trimmedInfo(iElem).trim(2);
    idNurbs = trimmedInfo(iElem).idNurbs;
    aNurbs = nurbs(idNurbs);
    vertCoord = Xe(1:3,:);

    % Quadrature for a NEFEM reference element
    [gauss,weights] = nefemQuad2DElementLocalCoordinates...
        (vertCoord,aNurbs,u1,u2,nIPNurbsEdge,nIPinterior);
    ngauss = length(weights);

    % Vandermonde matrix (nodes adapted)
%     nodesAdapted = nefemInterp2DAdaptedNodesElement(coordRef,vertCoord,aNurbs,u1,u2);
    nodesAdapted = inverseLinearMapping(vertCoord,Xe);
    V = Vandermonde_LP(nDeg,nodesAdapted);
    invV = inv(V');
    xe = Xe(:,1);
    ye = Xe(:,2);

    % Jacobian (computed with FEM reference element: straight faces)
%         v1 = Xe(1,:);  v2 = Xe(2,:);  v3 = Xe(3,:);
%         Jac = [(v2-v1)/2 ; (v3-v1)/2];
%         invJ = inv(Jac);
%         detJ = det(Jac);
%         dvolu = weights*detJ;

    % Loop in Gauss points
    nOfElementNodes = size(Xe,1);
    Ke = zeros(nOfElementNodes,nOfElementNodes);
    fe = zeros(nOfElementNodes,1);
    for g = 1:ngauss

        % Shape functions and derivatives (reference element)
        [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(gauss(g,:),nDeg);
        N_g = (invV*p)';
        Nxi_g = (invV*p_xi)';
        Neta_g = (invV*p_eta)';

        % Jacobian (computed with adaptedNodes)
        Jac = [Nxi_g*xe	Nxi_g*ye
            Neta_g*xe  Neta_g*ye];
        detJ = det(Jac);
        invJ = inv(Jac);
        dvolu_g = weights(g)*detJ;

        % Shape functions derivatives (cartesian element)
        Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
        Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;

        %phase and group celerities at current gauss point
        ccg_g = N_g*ccge;
        ccgx = Nx_g*ccge;
        ccgy = Ny_g*ccge;
        ccgGrad = [ccgx ccgy];
        k_g = N_g*ke;

        %Incident potential and derivative at the current gauss point
        xy_g = N_g*Xe;
        phi0_g = amp0*exp(sqrt(-1)*kvector*xy_g');
        phi0der_g = sqrt(-1)*kvector'*phi0_g;
        phi0Lap_g = -kvector*kvector'*phi0_g;

        %PML parameters
        sigma_g = N_g*sigmae;
        param = 1 + sqrt(-1)*sigma_g/omega;

        %Contribution of the current integration point to the elemental matrix and vector
%             dvolu_g = dvolu(g);
        K_e = ccg_g*((param(2)/param(1))*(Nx_g')*Nx_g + (param(1)/param(2))*(Ny_g')*Ny_g)*dvolu_g;
        M_e = ccg_g*k_g^2*param(1)*param(2)*(N_g')*N_g*dvolu_g;
        Ke = Ke + M_e - K_e;
        fe = fe + N_g'*param(1)*param(2)*(ccgGrad*phi0der_g + ccg_g*phi0Lap_g + ...
            k_g^2*ccg_g*phi0_g)*dvolu_g;
    end

    % Assembling
    f(Te) = f(Te) + fe;
    Te_transp = transpose(Te);
    aux_row = Te_transp(:,aux_ones);
    aux_col = Te(aux_ones,:);
    indK = (iElem-1)*nOfElementNodes^2+1:iElem*nOfElementNodes^2;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Kreal(indK) = real(Ke);
    Kimag(indK) = imag(Ke);
end

K = sparse(I,J,Kreal + sqrt(-1)*Kimag,nOfNodes,nOfNodes);
