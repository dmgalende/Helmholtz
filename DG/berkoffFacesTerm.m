function F = berkoffFacesTerm(U,mesh,theReferenceElement,infoFaces,...
    elementalMatricesInfo,t,incidentWave,BC,c,ccg)

nameNodal = mesh.fieldNames{mesh.indexElemPosCon(2)};
nameCon = mesh.fieldNames{mesh.indexElemPosCon(3)};
X = mesh.(nameNodal);
T = mesh.(nameCon);
nOfNodes = size(U,1);
nOfComponents = size(U,2);
elemFaceNodes = theReferenceElement.faceNodes;
nOfBoundaries = numel(mesh.boundaryNames);
N = theReferenceElement.N1d;
transN = transpose(N);
%%
faceMassReference = elementalMatricesInfo.faceMassReference;
ExteriorFaces = elementalMatricesInfo.ExteriorFaces;

% loop in Interior Faces
F = berkoffFacesTerm_Interior(X,T,U,theReferenceElement,infoFaces,faceMassReference,ccg);

% loop in Exterior Faces
for iboundary = 1:nOfBoundaries
    iname = mesh.boundaryNames{iboundary};
    icond = BC.values(iboundary);
    if icond == 1
        % NON Faces
        infoFaces_NON = infoFaces.(['exteriorFaces_' iname]);
        F = berkoff_NON_Faces(F,infoFaces_NON);
    elseif icond == 2
        % PEC condition
        iparam = BC.parameters{iboundary}{icond};
        infoFaces_PEC = infoFaces.(['exteriorFaces_' iname]);
        F = berkoff_PEC_Faces(F,infoFaces_PEC,iparam);
    elseif icond == 3
        % ABC condition
        iparam = BC.parameters{iboundary}{icond};
        infoFaces_ABC = infoFaces.(['exteriorFaces_' iname]);
        F = berkoff_ABC_Faces_FreeSpace(F,infoFaces_ABC,iparam);
    elseif icond == 4
        %Straight radiation boundary
        infoFaces_ABC = infoFaces.(['exteriorFaces_' iname]);
        F = berkoff_ABC_Faces_PML(F,infoFaces_ABC);
    end
end

    function F = berkoff_NON_Faces(F,infoFaces)
        %% Exterior faces: SYMMETRY
        nOfExteriorNON = size(infoFaces,1);
        for i = 1:nOfExteriorNON
            infoFace = infoFaces(i,:);
            iElemL = infoFace(1);
            % bottom informations
            nodesElem = T(iElemL,:); % element nodes in global numbering
            faceLocalNodes = elemFaceNodes(infoFace(2),:); % face nodes in local numbering Left Elem
            faceNodes = nodesElem(faceLocalNodes); % face nodes in global numbering Left Elem
            ccg_face = ccg(faceNodes);
            % Quadrature (weights) and outward unit normal (left)
            weightsFace = ExteriorFaces.(iname)(i).weights;
            n = ExteriorFaces.(iname)(i).extNormal;
            % Left and right states
            uLF = U(faceLocalNodes,:,iElemL);
            % Interpolate solution at gauss points
            uL = N*uLF;
            ccg_gauss = N*ccg_face;
            % Numerical minus normal flux-------------------
            FluxU_NON = berkoffFluxU_NON(uL,n,ccg_gauss);
            % Left element----------------------------------
            Fe = zeros(nOfNodes,nOfComponents);
            Ff = transN*(FluxU_NON.*repmat(weightsFace,1,nOfComponents));
            Fe(faceLocalNodes,:) = Ff;
            F(:,:,iElemL) = F(:,:,iElemL) + Fe;
        end
    end

    function F = berkoff_PEC_Faces(F,infoFaces,ABS)
        %% Exterior faces: PARTIAL ABSORPTION
        nOfExteriorPEC = size(infoFaces,1);
        for i = 1:nOfExteriorPEC
            infoFace = infoFaces(i,:);
            iElemL = infoFace(1);
            % bottom informations
            nodesElem = T(iElemL,:); % element nodes in global numbering
            faceLocalNodes = elemFaceNodes(infoFace(2),:); % face nodes in local numbering Left Elem
            faceNodes = nodesElem(faceLocalNodes); % face nodes in global numbering Left Elem
            ccg_face = ccg(faceNodes);
            c_face = c(faceNodes);
            % Quadrature (weights) and outward unit normal (left)
            weightsFace = ExteriorFaces.(iname)(i).weights;
            n = ExteriorFaces.(iname)(i).extNormal;
            gaussXYZ = ExteriorFaces.(iname)(i).gaussPoints;
            % Left and right states
            uLF = U(faceLocalNodes,:,iElemL);
            % Interpolate solution at gauss points
            uL = N*uLF;
            ccg_gauss = N*ccg_face;
            c_gauss = N*c_face;
            % Numerical minus normal flux-------------------
            FluxU_PEC = berkoffFluxU_PEC(uL,gaussXYZ,n,t,c_gauss,ccg_gauss,incidentWave,ABS);
            % Left element----------------------------------
            Fe = zeros(nOfNodes,nOfComponents);
            Ff = transN*(FluxU_PEC.*repmat(weightsFace,1,nOfComponents));
            Fe(faceLocalNodes,:) = Ff;
            F(:,:,iElemL) = F(:,:,iElemL) + Fe;
        end
    end

    function F = berkoff_ABC_Faces_FreeSpace(F,infoFaces,R)
        %% Exterior faces: ABC in free-space
        nOfExteriorABC = size(infoFaces,1);
        for i = 1:nOfExteriorABC
            infoFace = infoFaces(i,:);
            iElemL = infoFace(1);
            % bottom informations
            nodesElem = T(iElemL,:); % element nodes in global numbering
            faceLocalNodes = elemFaceNodes(infoFace(2),:); % face nodes in local numbering Left Elem
            faceNodes = nodesElem(faceLocalNodes); % face nodes in global numbering Left Elem
            ccg_face = ccg(faceNodes);
            c_face = c(faceNodes);
            % Quadrature (weights) and outward unit normal (left)
            weightsFace = ExteriorFaces.(iname)(i).weights;
            n = ExteriorFaces.(iname)(i).extNormal;
            % Left and right states
            uLF = U(faceLocalNodes,:,iElemL);
            % Interpolate solution at gauss points
            uL = N*uLF;
            ccg_gauss = N*ccg_face;
            c_gauss = N*c_face;
            % Numerical minus normal flux-------------------
            FluxU_ABC = berkoffFluxU_ABC(uL,n,c_gauss,ccg_gauss,R);
            % Left element----------------------------------
            Fe = zeros(nOfNodes,nOfComponents);
            Ff = transN*(FluxU_ABC.*repmat(weightsFace,1,nOfComponents));
            Fe(faceLocalNodes,:) = Ff;
            F(:,:,iElemL) = F(:,:,iElemL) + Fe;
        end
    end

    function F = berkoff_ABC_Faces_PML(F,infoFaces)
        %% Planar faces: ABC in PML
        nOfExteriorABC_PML = size(infoFaces,1);
        for i = 1:nOfExteriorABC_PML
            infoFace = infoFaces(i,:);
            iElemL = infoFace(1);
            % bottom informations
            nodesElem = T(iElemL,:); % element nodes in global numbering
            faceLocalNodes = elemFaceNodes(infoFace(2),:); % face nodes in local numbering Left Elem
            faceNodes = nodesElem(faceLocalNodes); % face nodes in global numbering Left Elem
            ccg_face = ccg(faceNodes);
            c_face = c(faceNodes);
            % Quadrature (weights) and outward unit normal (left)
            weightsFace = ExteriorFaces.(iname)(i).weights;
            n = ExteriorFaces.(iname)(i).extNormal;
            % Left and right states
            uLF = U(faceLocalNodes,:,iElemL);
            % Interpolate solution at gauss points
            uL = N*uLF;
            ccg_gauss = N*ccg_face;
            c_gauss = N*c_face;
            % Numerical minus normal flux-------------------
            FluxU_ABC = berkoffPMLFluxU_ABC(uL,n,c_gauss,ccg_gauss);
            % Left element----------------------------------
            Fe = zeros(nOfNodes,nOfComponents);
            Ff = transN*(FluxU_ABC.*repmat(weightsFace,1,nOfComponents));
            Fe(faceLocalNodes,:) = Ff;
            F(:,:,iElemL) = F(:,:,iElemL) + Fe;
        end
    end

end


