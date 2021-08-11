function elementalVectors  = computeAllElementalVectors(X,T,referenceElement,c,cg,incidentWave)

nOfElements = size(T,1);
nOfNodes = size(T,2);

elementalVectors.f1 = zeros(nOfNodes,nOfElements);
elementalVectors.f2 = zeros(nOfNodes,nOfElements);

for iElem = 1:nOfElements
    Te = T(iElem,:);
    coordElem = X(Te,:);
    c_el = c(Te);
    cg_el = cg(Te);
    [f1,f2] = computeElementalVectors(referenceElement, coordElem, c_el,cg_el,incidentWave);
    elementalVectors.f1(:,iElem) = f1;
    elementalVectors.f2(:,iElem) = f2;
end

function [f1,f2] = computeElementalVectors(referenceElement,coordElem,c,cg,incidentWave)

gaussWeights = referenceElement.IPweights;
nOfGauss = length(gaussWeights);
% incident wave computations
A = incidentWave.amplitude;
k = incidentWave.waveNumberValue;
omega = 2*pi/incidentWave.period;
angle = incidentWave.direction*(pi/180);
kvector = incidentWave.waveNumberValue*[cos(angle) sin(angle)];
kx = kvector(1); ky = kvector(2);
%

N = referenceElement.N;
Nxi  = referenceElement.Nxi;
Neta = referenceElement.Neta;

nOfNodes = size(referenceElement.NodesCoord,1);

f1 = zeros(nOfNodes,1);
f2 = zeros(nOfNodes,1);

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

    %coordinate of the Gauss point
    xg = N_g*xe;
    yg = N_g*ye;
    % c and cg at current gauss point
    c_g = N_g*c;
    cg_g = N_g*cg;
    
    alpha = c_g*cg_g; % ccg
    gamma = cg_g/c_g; % cg/c

    %alpha derivatives at current gauss point
    dalpha_dx = Nx_g*(c.*cg);
    dalpha_dy = Ny_g*(c.*cg);

    %     % Vectors
    gradAlphaK = dalpha_dx*kx + dalpha_dy*ky;
    f1 = f1 + N_g' * weight * A * ((alpha * k^2 - gamma*omega^2 )* cos(kx*xg + ky*yg) +...
        gradAlphaK * sin(kx*xg + ky*yg));
    f2 = f2 + N_g' * weight *A * ((alpha * k^2 - gamma*omega^2 )* sin(kx*xg + ky*yg) -...
        gradAlphaK * cos(kx*xg + ky*yg));

end