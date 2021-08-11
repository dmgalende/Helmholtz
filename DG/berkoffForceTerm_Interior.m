function B = berkoffForceTerm_Interior(F,U,theMesh,...
    t,incidentWave,referenceInvM,f1,f2,PML)

X = theMesh.X;
T = theMesh.T;
omega = 2*pi/incidentWave.period;
% initialization
B = zeros(size(U));

% Interior elements

for iElem = theMesh.interiorElements_FreeSpace
    Ue = U(:,:,iElem);
    Te = T(iElem,:);
    fe = f1(:,iElem)*cos(omega*t) + f2(:,iElem)*sin(omega*t);
    coordElem = X(Te,:);
    vertCoord = coordElem(1:3,:);
    J = Jacobian(vertCoord);
    detJ = det(J);
     % Contributions   
    B(:,1,iElem) = referenceInvM*( F(:,1,iElem)  + fe)/detJ;
    B(:,2,iElem) = referenceInvM* F(:,2,iElem)/detJ ;
    B(:,3,iElem) = referenceInvM* F(:,3,iElem)/detJ ;
    B(:,4,iElem) = referenceInvM* F(:,4,iElem)/detJ - Ue(:,1);

end

for iElem = theMesh.interiorElements_PML
    Ue = U(:,:,iElem);
    Te = T(iElem,:);    
    fe = f1(:,iElem)*cos(omega*t) + f2(:,iElem)*sin(omega*t);
    sigma_el = PML(Te,:);
    coordElem = X(Te,:);
    vertCoord = coordElem(1:3,:);
    J = Jacobian(vertCoord);
    detJ = det(J);
    PMLSourceTerm = berkoffPMLSource(sigma_el,Ue);
    % Contributions
    B(:,1,iElem) = referenceInvM*( F(:,1,iElem)  + fe)/detJ;
    B(:,2,iElem) = referenceInvM* F(:,2,iElem)/detJ ;
    B(:,3,iElem) = referenceInvM* F(:,3,iElem)/detJ ;
    B(:,4,iElem) = referenceInvM* F(:,4,iElem)/detJ - Ue(:,1);
    B(:,5,iElem) = referenceInvM* F(:,5,iElem)/detJ ;

    B(:,:,iElem) = B(:,:,iElem) + PMLSourceTerm;

end