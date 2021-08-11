function C = berkhoffNRBCMatrix_NEFEM...
    (X,TRobin,referenceElement,nurbs,trimmedInfo,k,ccg,R,NNZ_nefem)

%
% C = berkhoffNRBCMatrix_NEFEM...
%    (X,TRobin,referenceElement,nurbs,trimmedInfo,k,ccg,R)
%

%Number of elements and number of mesh nodes
nOfElements = size(TRobin,1);
nOfNodes = size(X,1);
nDeg = referenceElement.degree;

%Memory allocation
nOfElementNodes = size(TRobin,2);
aux_ones = ones(1,nOfElementNodes);
allocation = nOfElementNodes^2*nOfElements;
I = zeros(allocation,1);
J = I;
Creal = I;
Cimag = I;

% Gauss quadrature options
global nIPNurbsEdge

% FEM reference element
coordRef = referenceElement.NodesCoord;

% Vandermonde matrix
% V = Vandermonde_LP(nDeg,coordRef);
% invV = inv(V');
 
%Loop in 1D boundary elements
for iElem = 1:nOfElements
    Te = TRobin(iElem,:);
    Xe = X(Te,:); 
    ccge = ccg(Te);
    ke = k(Te);
    
    % NEFEM information
    u1 = trimmedInfo(iElem).trim(1);
    u2 = trimmedInfo(iElem).trim(2);
    idNurbs = trimmedInfo(iElem).idNurbs;
    aNurbs = nurbs(idNurbs);
    vertCoord = Xe(1:3,:);

    % Quadrature for a nurbs boundary integral
    [gauss,weights] = nefemQuad2DFaceParametricSpace...
        (aNurbs,u1,u2,nIPNurbsEdge);
    ngauss = length(weights);
    
    % Vandermonde matrix (nodes adapted)
%     nodesAdapted = nefemInterp2DAdaptedNodesElement(coordRef,vertCoord,aNurbs,u1,u2);
    nodesAdapted = inverseLinearMapping(vertCoord,Xe);
    V = Vandermonde_LP(nDeg,nodesAdapted);
    invV = inv(V');
    
    % Loop in Gauss points
    nOfElementNodes = size(Xe,1);
    Ce = zeros(nOfElementNodes,nOfElementNodes);    
    for g = 1:ngauss
        
        %Shape functions and derivatives at the current integration point
        igauss = gauss(g);
        pt = nurbsCurvePoint(aNurbs,igauss);
        gauss_xy = pt(1:2);
        gauss_xiEta = inverseLinearMapping(vertCoord,gauss_xy);
        p = orthopoly2D(gauss_xiEta,nDeg);
        N_g = (invV*p)';
        
        %Integration weight
        Cu = nurbsCurveDerivPoint(aNurbs,igauss);
        jacobianNurbs = norm(Cu(1:2));
        dline = weights(g)*jacobianNurbs;
        
        %phase and group celerities at current gauss point
        ccg_g = N_g*ccge;
        k_g = N_g*ke;
        
        %parameter for the Sommerfield NRB condition (radial boundary)
        param = sqrt(-1)*k_g - 1/(2*R);
  
        %Contribution of the current integration point to the elemental matrix
        Ce = Ce + ccg_g*param*(N_g')*N_g*dline;
    end
    
    % Assembling
    Te_transp = transpose(Te);
    aux_row = Te_transp(:,aux_ones);
    aux_col = Te(aux_ones,:);
    indK = (iElem-1)*nOfElementNodes^2+1:iElem*nOfElementNodes^2;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Creal(indK) = real(Ce);
    Cimag(indK) = imag(Ce);
end

C = sparse(I,J,Creal + sqrt(-1)*Cimag,nOfNodes,nOfNodes);
