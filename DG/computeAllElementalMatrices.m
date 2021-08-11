function elementalMatrices  = computeAllElementalMatrices(X,T,referenceElement)

nOfElements = size(T,1);
nOfNodes = size(T,2);

elementalMatrices.Cx = zeros(nOfNodes,nOfNodes,nOfElements);
elementalMatrices.Cy = zeros(nOfNodes,nOfNodes,nOfElements);
elementalMatrices.invM   = zeros(nOfNodes,nOfNodes,nOfElements);

for iElem = 1:nOfElements
    Te = T(iElem,:);
    coordElem = X(Te,:);    
    [M, Cx, Cy]  = computeElementalMatrices(referenceElement,coordElem);    
    elementalMatrices.invM(:,:,iElem)   = inv(M);
    elementalMatrices.Cx(:,:,iElem) = Cx;    
    elementalMatrices.Cy(:,:,iElem) = Cy;
end

function [M, Cx, Cy] = computeElementalMatrices...
    (theReferenceElement,coordElem)

gaussWeights = theReferenceElement.IPweights;
nOfGauss = length(gaussWeights);

N = theReferenceElement.N;
Nxi  = theReferenceElement.Nxi;
Neta = theReferenceElement.Neta;

nOfNodes = size(theReferenceElement.NodesCoord,1);
M = zeros(nOfNodes);
Cx = zeros(nOfNodes);
Cy = zeros(nOfNodes);

xe = coordElem(:,1);
ye = coordElem(:,2);

%Loop on Gauss points
for iGauss=1:nOfGauss
    N_g = N(iGauss,:);
    Nxi_g = Nxi(iGauss,:);
    Neta_g = Neta(iGauss,:);
    
    % Jacobian of the isoparametric transformation
    J = [Nxi_g*xe	  Nxi_g*ye
         Neta_g*xe    Neta_g*ye];
    
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end

    % Shape functions derivatives (cartesian element)
    invJ = inv(J);  % No necessary to invert!! Solve a sytem!!
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;
    weight = gaussWeights(iGauss)*det(J);
      
    % Matrices
    M =  M  + (weight*N_g')*N_g;
    Cx = Cx + (weight*N_g')*Nx_g;
    Cy = Cy + (weight*N_g')*Ny_g;

end