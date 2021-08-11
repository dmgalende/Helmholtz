function f = berkhoffReflectingVector_NEFEM...
    (X,T,referenceElement,nurbs,trimmedInfo,kvector,amp0,k,ccg,alpha)

%
% f = berkhoffReflectingVector_NEFEM...
%    (X,T,referenceElement,nurbs,trimmedInfo,kvector,amp0,k,ccg,alpha)
%

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfNodes = size(X,1);

%Memory allocation
f = zeros(nOfNodes,1);

% Gauss quadrature options
global nIPNurbsEdge

% FEM reference element
coordRef = referenceElement.NodesCoord;

% Vandermonde matrix
% nDeg = referenceElement.degree;
% V = Vandermonde_LP(nDeg,coordRef);
% invV = inv(V');
 
%Loop in 1D boundary elements
for iElem = 1:nOfElements
    Te = T(iElem,:);
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
    nDeg = referenceElement.degree;
    V = Vandermonde_LP(nDeg,nodesAdapted);
    invV = inv(V');
    
    % Loop in Gauss points
    nOfElementNodes = size(Xe,1);
    fe = zeros(nOfElementNodes,1); 
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
        
        %Incident potential at the current gauss point
        phi0_g = amp0*exp(sqrt(-1)*kvector*gauss_xy');
        phi0der_g = sqrt(-1)*kvector'*phi0_g;
        
        %Outward unit normal at the gauss point
        n_g = [Cu(2) -Cu(1)]/jacobianNurbs;
        
        %phase and group celerities at current gauss point
        ccg_g = N_g*ccge;
        k_g = N_g*ke;
        
        %Contribution of the current integration point to the elemental matrix
        fe = fe + ccg_g*(N_g')*(n_g*phi0der_g - sqrt(-1)*k_g*alpha*phi0_g)*dline;
    end
    
    % Assembing
    f(Te) = f(Te) + fe;
    clear fe
end