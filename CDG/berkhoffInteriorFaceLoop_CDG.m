function H = berkhoffInteriorFaceLoop_CDG...
    (X,T,theReferenceElement,infoFaces,alpha,sigma,omega,R,h)

nOfElements = size(T,1);
nOfElementNodes = size(theReferenceElement.NodesCoord,1);
nOfFaceNodes = numel(theReferenceElement.faceNodes1d);
nOfInteriorFaces = size(infoFaces.interiorFaces,1);
elemFaceNodes = theReferenceElement.faceNodes;% all faces nodes in local numbering
allocation = 3*nOfFaceNodes*nOfElementNodes*nOfElements + nOfElements*nOfElementNodes^2;

% H = spalloc(numel(T),numel(T),allocation);
% Hreal = spalloc(numel(T),numel(T),allocation);
% Himag = spalloc(numel(T),numel(T),allocation); %Non optimized NNZ for imag part!

I = zeros(allocation,1);
J = zeros(allocation,1);
Kr = zeros(allocation,1);
Ki = zeros(allocation,1);
index = 0;
element_ind =zeros(nOfElements,1) ;

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
aux_ones = ones(1,nOfElementNodes);
COUNT = 0.1;
time = cputime;
for i = 1:nOfInteriorFaces
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
    indF = (iElemF-1)*nOfElementNodes+1:iElemF*nOfElementNodes;
    indS = (iElemS-1)*nOfElementNodes+1:iElemS*nOfElementNodes;
    
    % calculate n1 (=-n2)
    [n,A] = faceNormalArea(infoFace,theReferenceElement,X,T);
    M_iFace = M_reference_face * A/2; % mass matrix of the current face (/2 because the lenght of the reference face is 2)

    %___________________
    % Compact DG method
    % Consistent switch or Natural switch (comment/uncomment)
    %___________________

    % Consistent switch C12
    if faceNodes(1)<faceNodes(end)
        C12=1;
    else
        C12=-1;
    end
    % Natural switch
    %     if iElemF<iElemS
    %         C12=1;
    %     else
    %         C12=-1;
    %     end
    beta = 0;
    C11 = beta/h;

    % for the lifting functions
    R_ElemF = R{iElemF};
    R_ElemS = R{iElemS};

    %     first element
    [FA11, FA12, FA21, FA22] = lifting(1,theReferenceElement,n,infoFace,R_ElemF,M_iFace,C12);
    % second element
    [SA11, SA12, SA21, SA22] = lifting(2,theReferenceElement,-n,infoFace,R_ElemS,M_iFace,-C12);


    %     %___________________
    %     % Interior Penalty Method (comment/uncomment previous switch lines) - no switch
    %     %___________________
    %
    %     C12 = zeros(size(n));
    %     C11_0 = 100; % to be calibrated with the first mesh
    %     C11 = C11_0/h;
    %
    %     FA11 = zeros(nOfElementNodes,nOfElementNodes);
    %     FA12 = FA11;
    %     FA21 = FA11;
    %     FA22 = FA11;
    %     SA11 = FA11;
    %     SA12 = FA11;
    %     SA21 = FA11;
    %     SA22 = FA11;

    % Face contribution
    [XX,XY,YX,YY] = faceContribution(theReferenceElement,XeF,XeS,n,A,infoFace,alpha_f,C11,C12,sigma_f,omega,...
        Nx_face,Ny_face);


    %     % OLD ASSEMBLING
    %     % assembly: 1st element
    %     %contribution of the First Elem to First Elem
    %     H(indF,indF) = H(indF,indF) + XX + FA11 + SA22;
    %     %contribution of the Second Elem to First Elem
    %     H(indF,indS) = H(indF,indS) + XY + FA12 + SA21;
    %
    %     % assembly: 2nd element
    %     %contribution of the First Elem to Second Elem
    %     H(indS,indF) = H(indS,indF) + YX + FA21 + SA12;
    %     %contribution of the Second Elem to Second Elem
    %     H(indS,indS) = H(indS,indS) + YY + FA22 + SA11;
    %     % ASSEMBLING WITH SSO
    %     if ~isreal(R_ElemF) || ~isreal(R_ElemS)
    %         %assembly: 1st element
    %         %contribution of the First Elem to First Elem
    %         SSO( Himag, imag(XX + FA11 + SA22), indF, indF, '+' )
    %         SSO( Hreal, real(XX + FA11 + SA22), indF, indF, '+' )
    %         %contribution of the Second Elem to First Elem
    %         SSO( Himag, imag(XY + FA12 + SA21), indF, indS, '+' )
    %         SSO( Hreal, real(XY + FA12 + SA21), indF, indS, '+' )
    %
    %         %assembly: 2nd element
    %         %contribution of the First Elem to Second Elem
    %         SSO( Himag, imag(YX + FA21 + SA12), indS, indF, '+' )
    %         SSO( Hreal, real(YX + FA21 + SA12), indS, indF, '+' )
    %         %contribution of the Second Elem to Second Elem
    %         SSO( Himag, imag(YY + FA22 + SA11), indS, indS, '+' )
    %         SSO( Hreal, real(YY + FA22 + SA11), indS, indS, '+' )
    %     else
    %         %assembly: 1st element
    %         %contribution of the First Elem to First Elem
    %         SSO( Hreal, XX + FA11 + SA22, indF, indF, '+' )
    %         %contribution of the Second Elem to First Elem
    %         SSO( Hreal, XY + FA12 + SA21, indF, indS, '+' )
    %
    %         %assembly: 2nd element
    %         %contribution of the First Elem to Second Elem
    %         SSO( Hreal, YX + FA21 + SA12, indS, indF, '+' )
    %         %contribution of the Second Elem to Second Elem
    %         SSO( Hreal, YY + FA22 + SA11, indS, indS, '+' )
    %     end

    % NEW ASSEMBLING
    indF_transp = transpose(indF);
    indS_transp = transpose(indS);
    aux_mat = zeros(nOfElementNodes^2,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indF_transp(:,aux_ones);
    aux_col = indF(aux_ones,:);
    aux_mat(:) =  XX + FA11 + SA22;
    if element_ind(iElemF)==0
        index = index(end) + (1:nOfElementNodes^2);
        I(index) = aux_row(:);
        J(index) = aux_col(:);
        Kr(index) = real(aux_mat);
        Ki(index) = imag(aux_mat);
        element_ind(iElemF) = index(1);
    else
        prev_ind = element_ind(iElemF)-1+ (1:nOfElementNodes^2);
        Kr(prev_ind) = Kr(prev_ind) + real(aux_mat);
        Ki(prev_ind) = Ki(prev_ind) + imag(aux_mat);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indF_transp(:,aux_ones); aux_row = aux_row(:);
    aux_col = indS(aux_ones,:); aux_col = aux_col(:);
    aux_mat(:) =   XY + FA12 + SA21;
    aux_ind = aux_mat~=0;
    index = index(end) + (1:length(aux_mat(aux_ind)));
    I(index) = aux_row(aux_ind);
    J(index) = aux_col(aux_ind);
    Kr(index) =  real(aux_mat(aux_ind));
    Ki(index) =  imag(aux_mat(aux_ind));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indS_transp(:,aux_ones); aux_row = aux_row(:);
    aux_col = indF(aux_ones,:); aux_col = aux_col(:);
    aux_mat(:) = YX + FA21 + SA12;
    aux_ind = aux_mat~=0;
    index = index(end) + (1:length(aux_mat(aux_ind)));
    I(index) = aux_row(aux_ind);
    J(index) = aux_col(aux_ind);
    Kr(index) =  real(aux_mat(aux_ind));
    Ki(index) =  imag(aux_mat(aux_ind));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    aux_row = indS_transp(:,aux_ones);
    aux_col = indS(aux_ones,:);
    aux_mat(:) = YY + FA22 + SA11;
    if element_ind(iElemS)==0
        index = index(end) + (1:nOfElementNodes^2);
        I(index) = aux_row(:);
        J(index) = aux_col(:);
        Kr(index) = real(aux_mat);
        Ki(index) = imag(aux_mat);
        element_ind(iElemS) = index(1);
    else
        prev_ind = element_ind(iElemS)-1+ (1:nOfElementNodes^2);
        Kr(prev_ind) = Kr(prev_ind) + real(aux_mat);
        Ki(prev_ind) = Ki(prev_ind) + imag(aux_mat);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i > nOfInteriorFaces*COUNT
        disp(['Calculated: ' num2str(100*COUNT) '% ; time elapsed ' num2str(cputime-time)])
        time = cputime;
        COUNT = COUNT + 0.1;
    end
end
% FOR THE SSO ASSEMBLING
% H = Hreal + sqrt(-1)*Himag;
disp('Creating sparse matrix')
disp(' ')
H = sparse(I(I~=0),J(I~=0),Kr(I~=0)+sqrt(-1)*Ki(I~=0),numel(T),numel(T));

