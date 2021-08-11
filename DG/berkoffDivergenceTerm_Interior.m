function divF = berkoffDivergenceTerm_Interior(F,U,...
    elementalMatricesReference,theMesh,ccg)

X = theMesh.X;
T = theMesh.T;
% initialization
divF = zeros(size(U));

% Planar elements in free space
Cxi = elementalMatricesReference.Cxi;
Ceta = elementalMatricesReference.Ceta;

for iElem = theMesh.interiorElements_FreeSpace
    Ue = U(:,:,iElem);
    Te = T(iElem,:);
    vertCoord = X(Te(1:3),:);
    J = Jacobian(vertCoord);
    invJ = inv(J);
    detJ = det(J);
    Cx = invJ(1,1)*Cxi + invJ(1,2)*Ceta;
    Cy = invJ(2,1)*Cxi + invJ(2,2)*Ceta;
    Cx_transp =  Cx';
    Cy_transp =  Cy';
    ccg_el = ccg(Te);
    % Contributions
    divF(:,1,iElem) =  F(:,1,iElem) + detJ*(Cx_transp*(ccg_el.*Ue(:,2)) + Cy_transp*(ccg_el.*Ue(:,3)));
    divF(:,2,iElem) = F(:,2,iElem) + detJ*(Cx_transp* Ue(:,1)) ;
    divF(:,3,iElem) =  F(:,3,iElem) +detJ*(Cy_transp * Ue(:,1));
    divF(:,4,iElem) =  F(:,4,iElem);

end

for iElem = theMesh.interiorElements_PML
    Ue = U(:,:,iElem);
    Te = T(iElem,:);
    vertCoord = X(Te(1:3),:);
    J = Jacobian(vertCoord);
    invJ = inv(J);
    detJ = det(J);
    Cx = invJ(1,1)*Cxi + invJ(1,2)*Ceta;
    Cy = invJ(2,1)*Cxi + invJ(2,2)*Ceta;
    Cx_transp =  Cx';
    Cy_transp =  Cy';
    ccg_el = ccg(Te);
    % Contributions
    divF(:,1,iElem) =  F(:,1,iElem) + detJ*(Cx_transp*(ccg_el.*Ue(:,2)) + Cy_transp*(ccg_el.*Ue(:,3)));
    divF(:,2,iElem) = F(:,2,iElem) + detJ*(Cx_transp* Ue(:,1)) ;
    divF(:,3,iElem) =  F(:,3,iElem) +detJ*(Cy_transp * Ue(:,1));
    divF(:,4,iElem) =  F(:,4,iElem);
    %PML components
    divF(:,5,iElem) =  F(:,5,iElem) + detJ*Cy_transp*(ccg_el.*Ue(:,3));

end
