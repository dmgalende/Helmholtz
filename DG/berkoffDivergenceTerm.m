function divF = berkoffDivergenceTerm(F,U,theMesh,...
    elementalMatricesInfo,ccg)

T = theMesh.T;
nOfInteriorElements = numel(theMesh.interiorElements_FreeSpace) + ...
                      numel(theMesh.interiorElements_PML);
%--------------------------------------------------------- Interior Elements
divF = berkoffDivergenceTerm_Interior(F,U,...
     elementalMatricesInfo.elementalMatricesReference,theMesh,ccg);

%---------------------------------------------------------- Exterior elements

Cx = elementalMatricesInfo.elementalMatricesExterior.Cx;
Cy = elementalMatricesInfo.elementalMatricesExterior.Cy;


for iElem = theMesh.exteriorElements_FreeSpace
    Ue = U(:,:,iElem);
    ccg_el = ccg(T(iElem,:));
    jElem = iElem - nOfInteriorElements;
    Cx_Elem_transp =  Cx(:,:,jElem)';
    Cy_Elem_transp =  Cy(:,:,jElem)';
    
    divF(:,1,iElem) =  F(:,1,iElem) + Cx_Elem_transp*(ccg_el.*Ue(:,2)) + ...
                                      Cy_Elem_transp*(ccg_el.*Ue(:,3)); 
    divF(:,2,iElem) = F(:,2,iElem) + Cx_Elem_transp* Ue(:,1) ; 
    divF(:,3,iElem) =  F(:,3,iElem) +Cy_Elem_transp * Ue(:,1); 
    divF(:,4,iElem) =  F(:,4,iElem);

end

for iElem = theMesh.exteriorElements_PML
    Ue = U(:,:,iElem);
    ccg_el = ccg(T(iElem,:));
    jElem = iElem - nOfInteriorElements;    
    Cx_Elem_transp =  Cx(:,:,jElem)';
    Cy_Elem_transp =  Cy(:,:,jElem)';

    divF(:,1,iElem) =  F(:,1,iElem) + Cx_Elem_transp*(ccg_el.*Ue(:,2)) + ...
        Cy_Elem_transp*(ccg_el.*Ue(:,3));
    divF(:,2,iElem) = F(:,2,iElem) + Cx_Elem_transp* Ue(:,1) ;
    divF(:,3,iElem) =  F(:,3,iElem) +Cy_Elem_transp * Ue(:,1);
    divF(:,4,iElem) =  F(:,4,iElem);
    divF(:,5,iElem) =  F(:,5,iElem) + Cy_Elem_transp*(ccg_el.*Ue(:,3));    %
end