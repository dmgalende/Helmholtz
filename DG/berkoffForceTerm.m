function B = berkoffForceTerm(F,U,theMesh,...
    elementalMatricesInfo,t,incidentWave,PML)

T = theMesh.T;
invM = elementalMatricesInfo.elementalMatricesExterior.invM;
omega = 2*pi/incidentWave.period;
% modulation
aux = omega*t/(2*pi);
if aux<=1
    modulation = aux;
else
    modulation = 1;
end;
%

f1 = modulation * elementalMatricesInfo.elementalVectors.f1;
f2 = modulation * elementalMatricesInfo.elementalVectors.f2;
nOfInteriorElements = numel(theMesh.interiorElements_FreeSpace) + ...
                      numel(theMesh.interiorElements_PML);
referenceInvM = elementalMatricesInfo.elementalMatricesReference.invM;

%----------------------------- Interior Elements
B = berkoffForceTerm_Interior(F,U,theMesh,...
    t,incidentWave,referenceInvM,f1,f2,PML);

%----------------------------- Exterior Elements

for iElem = theMesh.exteriorElements_FreeSpace
    Ue = U(:,:,iElem);
    fe = f1(:,iElem)*cos(omega*t) + f2(:,iElem)*sin(omega*t);
    % Contributions
    jElem = iElem - nOfInteriorElements;
    invM_Elem = invM(:,:,jElem);
    
    B(:,1,iElem) = invM_Elem*( F(:,1,iElem)  + fe) ;
    B(:,2,iElem) = invM_Elem* F(:,2,iElem) ;
    B(:,3,iElem) = invM_Elem* F(:,3,iElem) ;
    B(:,4,iElem) = invM_Elem* F(:,4,iElem) - Ue(:,1);
    
end

for iElem = theMesh.exteriorElements_PML
    Ue = U(:,:,iElem);
    Te = T(iElem,:);
    sigma_el = PML(Te,:);
    fe = f1(:,iElem)*cos(omega*t) + f2(:,iElem)*sin(omega*t);
    % Contributions
    jElem = iElem - nOfInteriorElements;
    invM_Elem = invM(:,:,jElem);
    PMLSourceTerm = berkoffPMLSource(sigma_el,Ue);
    
    B(:,1,iElem) = invM_Elem*( F(:,1,iElem)  + fe);
    B(:,2,iElem) = invM_Elem* F(:,2,iElem) ;
    B(:,3,iElem) = invM_Elem* F(:,3,iElem) ;
    B(:,4,iElem) = invM_Elem* F(:,4,iElem) - Ue(:,1);
    B(:,5,iElem) = invM_Elem* F(:,5,iElem) ;

    B(:,:,iElem) = B(:,:,iElem) + PMLSourceTerm;


end