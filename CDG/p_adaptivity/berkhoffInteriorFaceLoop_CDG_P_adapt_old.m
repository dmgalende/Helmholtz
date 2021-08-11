function H = berkhoffInteriorFaceLoop_CDG_P_adapt_old...
    (X,T,theReferenceElement,infoFaces,alpha,sigma,omega,R,p_adaptivity)

nOfElements = size(T,1);
elem_p = p_adaptivity.elem_p;
elem_check = p_adaptivity.elem_check;
nOfElementNodes = size(theReferenceElement.NodesCoord,1);
nOfElementNodes_vect = 0.5*((elem_p+1).*(elem_p+2));
nOfFaceNodes_vect = elem_p+1;
nOfInteriorFaces = size(infoFaces.interiorFaces,1);
interiorFaces = 1:nOfInteriorFaces;
elemFaceNodes = theReferenceElement.faceNodes;% all faces nodes in local numbering

% Allocation
allocation = 3*sum(nOfFaceNodes_vect.*nOfElementNodes_vect) + sum(nOfElementNodes_vect.^2); 
dim = 0.5*sum((elem_p+1).*(elem_p+2));
I = zeros(allocation,1);
J = zeros(allocation,1);
Kr = zeros(allocation,1);
Ki = zeros(allocation,1);

% indexes for the assembly process
index = 0;
element_ind =zeros(nOfElements,1);

% faces involving refined elements
face_check = elem_check(infoFaces.interiorFaces(:,1)) | elem_check(infoFaces.interiorFaces(:,3));

% face mass matrix and shape function for faces not belonging to refined
% elements
M_reference_face = faceMassMatrix(theReferenceElement); % mass matrix of the reference face
Nx_face = zeros(nOfElementNodes,numel(theReferenceElement.IPcoordinates1d),3);
Ny_face = Nx_face ;
for iFace = 1:3
    coordGaussPointsFaceRefEl = one2twoDmapping(iFace, theReferenceElement);
    [shapeFun, shapeFunDer] = evalContinuousShapeFunctionsAtTriPoints(coordGaussPointsFaceRefEl,...
        theReferenceElement.degree, theReferenceElement.NodesCoord);
    Nx_face(:,:,iFace) = shapeFunDer(:,:,1);
    Ny_face(:,:,iFace) = shapeFunDer(:,:,2);
end

%% NON-Refined Elements
for i = interiorFaces(~face_check)
    infoFace = infoFaces.interiorFaces(i,:);
    iElemF = infoFace(1);
    iElemS = infoFace(3);
    
    % elements coordinates
    XeF = X(T(iElemF,:),:);
    XeS = X(T(iElemS,:),:);
    
    %ccg coefficient in the face nodes
    nodesElem = T(iElemF,:); % element nodes in global numbering
    faceNodes = nodesElem(elemFaceNodes(infoFace(2),:)); % face nodes in global numbering
    alpha_f = alpha(faceNodes);
    sigma_f = sigma(faceNodes,:);
    
    % nodal indices for the first and second element
    indF = 0.5*sum((elem_p(1:iElemF-1)+1).*(elem_p(1:iElemF-1)+2))+1 : ...
        0.5*sum((elem_p(1:iElemF)+1).*(elem_p(1:iElemF)+2));  % index for the assembly process
    indS = 0.5*sum((elem_p(1:iElemS-1)+1).*(elem_p(1:iElemS-1)+2))+1 : ...
        0.5*sum((elem_p(1:iElemS)+1).*(elem_p(1:iElemS)+2));  % index for the assembly process
    
    % calculate n1 (=-n2)
    [n,A] = faceNormalArea(infoFace,theReferenceElement,X,T);
    M_iFace = M_reference_face * A/2; % mass matrix of the current face (/2 because the lenght of the reference face is 2)
    
    % Consistent switch C12
    if faceNodes(1)<faceNodes(end)
        C12=1;
    else
        C12=-1;
    end
    C11 = 0;
    
    % for the lifting functions
    R_ElemF = R{iElemF};
    R_ElemS = R{iElemS};
    
    % first element
    [FA11, FA12, FA21, FA22] = lifting(1,theReferenceElement,n,infoFace,R_ElemF,M_iFace,C12);
    
    % second element
    [SA11, SA12, SA21, SA22] = lifting(2,theReferenceElement,-n,infoFace,R_ElemS,M_iFace,-C12);

    % Face contribution
    [XX,XY,YX,YY] = faceContribution(theReferenceElement,XeF,XeS,n,A,infoFace,alpha_f,C11,C12,sigma_f,omega,...
        Nx_face,Ny_face);

    % NEW ASSEMBLING
    indF_transp = transpose(indF);
    indS_transp = transpose(indS);
    aux_mat = zeros(nOfElementNodes^2,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indF_transp(:,ones(1,nOfElementNodes)); aux_row = aux_row(:);
    aux_col = indF(ones(1,nOfElementNodes),:); aux_col = aux_col(:);
    aux_mat(:) =  XX + FA11 + SA22;
    if element_ind(iElemF)==0
        index = index(end) + (1:nOfElementNodes^2);
        I(index) = aux_row;
        J(index) = aux_col;
        Kr(index) = real(aux_mat);
        Ki(index) = imag(aux_mat);
        element_ind(iElemF) = index(1);
    else
        prev_ind = element_ind(iElemF)-1+ (1:nOfElementNodes^2);
        Kr(prev_ind) = Kr(prev_ind) + real(aux_mat);
        Ki(prev_ind) = Ki(prev_ind) + imag(aux_mat);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indF_transp(:,ones(1,nOfElementNodes)); aux_row = aux_row(:);
    aux_col = indS(ones(1,nOfElementNodes),:); aux_col = aux_col(:);
    aux_mat(:) =   XY + FA12 + SA21;
    aux_ind = aux_mat~=0;
    index = index(end) + (1:length(aux_mat(aux_ind)));
    I(index) = aux_row(aux_ind);
    J(index) = aux_col(aux_ind);
    Kr(index) =  real(aux_mat(aux_ind));
    Ki(index) =  imag(aux_mat(aux_ind));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indS_transp(:,ones(1,nOfElementNodes)); aux_row = aux_row(:);
    aux_col = indF(ones(1,nOfElementNodes),:); aux_col = aux_col(:);
    aux_mat(:) = YX + FA21 + SA12;
    aux_ind = aux_mat~=0;
    index = index(end) + (1:length(aux_mat(aux_ind)));
    I(index) = aux_row(aux_ind);
    J(index) = aux_col(aux_ind);
    Kr(index) =  real(aux_mat(aux_ind));
    Ki(index) =  imag(aux_mat(aux_ind));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indS_transp(:,ones(1,nOfElementNodes)); aux_row = aux_row(:);
    aux_col = indS(ones(1,nOfElementNodes),:); aux_col = aux_col(:);
    aux_mat(:) = YY + FA22 + SA11;
    if element_ind(iElemS)==0
        index = index(end) + (1:nOfElementNodes^2);
        I(index) = aux_row;
        J(index) = aux_col;
        Kr(index) = real(aux_mat);
        Ki(index) = imag(aux_mat);
        element_ind(iElemS) = index(1);
    else
        prev_ind = element_ind(iElemS)-1+ (1:nOfElementNodes^2);
        Kr(prev_ind) = Kr(prev_ind) + real(aux_mat);
        Ki(prev_ind) = Ki(prev_ind) + imag(aux_mat);
    end
end

%% Refined Elements
for i = interiorFaces(face_check)
    infoFace = infoFaces.interiorFaces(i,:);
    if elem_p(infoFace(1)) >= elem_p(infoFace(3));
        iElemF = infoFace(1);
        iFaceF = infoFace(2);
        iElemS = infoFace(3);
        iFaceS = infoFace(4);
    else
        iElemF = infoFace(3);
        iFaceF = infoFace(4);
        iElemS = infoFace(1);
        iFaceS = infoFace(2);
        infoFace = infoFace([3 4 1 2 5]);
    end
    
    % first (refined) element
    p_F = elem_p(iElemF);
    refEl_F = p_adaptivity.referenceElements.(['P' num2str(p_F)]);
    nOfElementNodes_F = 0.5*(p_F+1)*(p_F+2);
    indRef_F = 0.5*sum((elem_p(elem_check(1:iElemF-1))+1).*(elem_p(elem_check(1:iElemF-1))+2))+(1:nOfElementNodes_F);
    Xe_F = p_adaptivity.coord(indRef_F,:);
    faceNodes = refEl_F.faceNodes;

    % second (refined or not) element
    p_S = elem_p(iElemS);
    nOfElementNodes_S = 0.5*(p_S+1)*(p_S+2);
    if isfield(p_adaptivity.referenceElements,['P' num2str(p_S)]);
        % refined
        refEl_S = p_adaptivity.referenceElements.(['P' num2str(p_S)]);
        indRef_S = 0.5*sum((elem_p(elem_check(1:iElemS-1))+1).*(elem_p(elem_check(1:iElemS-1))+2))+(1:nOfElementNodes_S);
        Xe_S = p_adaptivity.coord(indRef_S,:);
    else
        % not refined
        refEl_S = theReferenceElement;
        Xe_S = X(T(iElemS,:),:);
    end
    
    %ccg coefficient in the face nodes for the highest P
    alpha_el = p_adaptivity.ccg(indRef_F);
    sigma_el = p_adaptivity.sigma(indRef_F,:); 
    alpha_f = alpha_el(refEl_F.faceNodes(iFaceF,:));
    sigma_f = sigma_el(refEl_F.faceNodes(iFaceF,:),:);
    
    % nodal indices for the first and second element
    indF = 0.5*sum((elem_p(1:iElemF-1)+1).*(elem_p(1:iElemF-1)+2))+1 : ...
        0.5*sum((elem_p(1:iElemF)+1).*(elem_p(1:iElemF)+2));  % index for the assembly process
    indS = 0.5*sum((elem_p(1:iElemS-1)+1).*(elem_p(1:iElemS-1)+2))+1 : ...
        0.5*sum((elem_p(1:iElemS)+1).*(elem_p(1:iElemS)+2));  % index for the assembly process
    
    % face mass matrix and shape function for faces belonging to refined
    % elements
    M_reference_face_F = faceMassMatrix(refEl_F); % mass matrix of the reference face
    M_reference_face_S = faceMassMatrix(refEl_S); % mass matrix of the reference face
    coordGaussPointsFaceRefEl_F = one2twoDmapping(iFaceF, refEl_F);
    coordGaussPointsFaceRefEl_S = one2twoDmapping(iFaceS, refEl_F); % I use the same refEl to have the same gauss points
    [shapeFun_F, shapeFunDer_F] = evalContinuousShapeFunctionsAtTriPoints(coordGaussPointsFaceRefEl_F,...
        refEl_F.degree, refEl_F.NodesCoord);
    [shapeFun_S, shapeFunDer_S] = evalContinuousShapeFunctionsAtTriPoints(coordGaussPointsFaceRefEl_S,...
        refEl_S.degree, refEl_S.NodesCoord);
    N_face_F = shapeFun_F'; N_face_F = N_face_F(:,refEl_F.faceNodes(iFaceF,:));
    N_face_S = shapeFun_S'; N_face_S = N_face_S(:,refEl_S.faceNodes(iFaceS,:));
    Nx_face_F = shapeFunDer_F(:,:,1)';
    Ny_face_F = shapeFunDer_F(:,:,2)';
    Nx_face_S = shapeFunDer_S(:,:,1)';
    Ny_face_S = shapeFunDer_S(:,:,2)';

    % calculate n1 (=-n2)
    [n,A] = faceNormalArea(infoFace,theReferenceElement,X,T);
    M_iFace_F = M_reference_face_F * A/2; % mass matrix of the current face (/2 because the lenght of the reference face is 2)
    M_iFace_S = M_reference_face_S * A/2;
    M_face_12 = A/2*faceMassMatrix_x(N_face_F,N_face_S,refEl_F,refEl_S);
    
    % Consistent switch C12
    if faceNodes(1)<faceNodes(end)
        C12=1;
    else
        C12=-1;
    end
    C11 = 0;
    
    % for the lifting functions
    R_ElemF = R{iElemF};
    R_ElemS = R{iElemS};
    
    % first element
    [FA11, FA12, FA21, FA22] = lifting_P_adapt(iFaceF,iFaceS,refEl_F,refEl_S,n,R_ElemF,M_iFace_F,M_face_12,C12);
    
    % second element
    [SA11, SA12, SA21, SA22] = lifting_P_adapt(iFaceS,iFaceF,refEl_S,refEl_F,-n,R_ElemS,M_iFace_S,M_face_12',-C12);

    % Face contribution
    [XX,XY,YX,YY] = faceContribution_P_adapt(iFaceF,iFaceS,refEl_F,refEl_S,Xe_F,Xe_S,n,A,alpha_f,...
        C11,C12,sigma_f,omega,N_face_F,N_face_S,Nx_face_F,Nx_face_S,Ny_face_F,Ny_face_S);

    % NEW ASSEMBLING
    indF_transp = transpose(indF);
    indS_transp = transpose(indS);
    aux_mat_1 = zeros(nOfElementNodes_F^2,1);
    aux_mat_12 = zeros(nOfElementNodes_F*nOfElementNodes_S,1);
    aux_mat_2 = zeros(nOfElementNodes_S^2,1);
    aux_ones_1 = ones(1,nOfElementNodes_F);
    aux_ones_2 = ones(1,nOfElementNodes_S);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indF_transp(:,aux_ones_1); aux_row = aux_row(:);
    aux_col = indF(aux_ones_1,:); aux_col = aux_col(:);
    aux_mat_1(:) =  XX + FA11 + SA22;
    if element_ind(iElemF)==0
        index = index(end) + (1:nOfElementNodes_F^2);
        I(index) = aux_row;
        J(index) = aux_col;
        Kr(index) = real(aux_mat_1);
        Ki(index) = imag(aux_mat_1);
        element_ind(iElemF) = index(1);
    else
        prev_ind = element_ind(iElemF)-1+ (1:nOfElementNodes_F^2);
        Kr(prev_ind) = Kr(prev_ind) + real(aux_mat_1);
        Ki(prev_ind) = Ki(prev_ind) + imag(aux_mat_1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indF_transp(:,aux_ones_2); aux_row = aux_row(:);
    aux_col = indS(aux_ones_1,:); aux_col = aux_col(:);
    aux_mat_12(:) =   XY + FA12 + SA21;
    aux_ind = aux_mat_12~=0;
    index = index(end) + (1:length(aux_mat_12(aux_ind)));
    I(index) = aux_row(aux_ind);
    J(index) = aux_col(aux_ind);
    Kr(index) =  real(aux_mat_12(aux_ind));
    Ki(index) =  imag(aux_mat_12(aux_ind));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indS_transp(:,aux_ones_1); aux_row = aux_row(:);
    aux_col = indF(aux_ones_2,:); aux_col = aux_col(:);
    aux_mat_12(:) = YX + FA21 + SA12;
    aux_ind = aux_mat_12~=0;
    index = index(end) + (1:length(aux_mat_12(aux_ind)));
    I(index) = aux_row(aux_ind);
    J(index) = aux_col(aux_ind);
    Kr(index) =  real(aux_mat_12(aux_ind));
    Ki(index) =  imag(aux_mat_12(aux_ind));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indS_transp(:,aux_ones_2); aux_row = aux_row(:);
    aux_col = indS(aux_ones_2,:); aux_col = aux_col(:);
    aux_mat_2(:) = YY + FA22 + SA11;
    if element_ind(iElemS)==0
        index = index(end) + (1:nOfElementNodes_S^2);
        I(index) = aux_row;
        J(index) = aux_col;
        Kr(index) = real(aux_mat_2);
        Ki(index) = imag(aux_mat_2);
        element_ind(iElemS) = index(1);
    else
        prev_ind = element_ind(iElemS)-1+ (1:nOfElementNodes_S^2);
        Kr(prev_ind) = Kr(prev_ind) + real(aux_mat_2);
        Ki(prev_ind) = Ki(prev_ind) + imag(aux_mat_2);
    end
end

% create sparse matrix
disp('Creating sparse matrix')
disp(' ')
H = sparse(I(I~=0),J(I~=0),Kr(I~=0)+sqrt(-1)*Ki(I~=0),dim,dim);
