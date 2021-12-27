function [Kint,Mint,Kxpml,Kypml,Mpml] = PGDberkhoffVolumeMatrices...
    (X,T,theReferenceElement,vcoef1,vcoef2,vcoef3,vcoef4,inPMLelems)

%Vectorized elements
nElemsVectorized = 1024;

%Number of elements, number of mesh nodes, initialization of coefficients
nOfElements = size(T,1);
nOfNodes = size(X,1);
nOfTerms1 = size(vcoef1,2);
nOfTerms2 = size(vcoef2,2);
nOfTerms3 = size(vcoef3,2);
nOfTerms4 = size(vcoef4,2);
coef1e = cell(nOfTerms1,1);
coef2e = cell(nOfTerms2,1);
coef3e = cell(nOfTerms3,1);
coef4e = cell(nOfTerms4,1);

%Initialize output cells
Kint  = coef1e;
Mint  = coef4e;
Kxpml = coef2e;
Kypml = coef3e;
Mpml  = coef4e;

%Initialize allocation cells
Kint_real  = coef1e;
Kint_imag  = coef1e;
Mint_real  = coef4e;
Mint_imag  = coef4e;
Kxpml_real = coef2e;
Kxpml_imag = coef2e;
Kypml_real = coef3e;
Kypml_imag = coef3e;
Mpml_real  = coef4e;
Mpml_imag  = coef4e;

%Information of the reference element
IPw = theReferenceElement.IPweights;                % ng??1
N = theReferenceElement.N;                          % ng??nOfElementNodes
Nxi = theReferenceElement.Nxi;                      % ng??nOfElementNodes
Neta = theReferenceElement.Neta;                    % ng??nOfElementNodes
ngauss = length(IPw);

%Memory allocation for the interior elements
nOfElementNodes = size(T,2);
aux_ones = ones(1,nOfElementNodes);
allocation = nOfElementNodes^2*nOfElements;
Iint = zeros(allocation,1);                            % allocation??1
Jint = Iint;                                              % allocation??1
Kint_real(:) = {Iint};
Kint_imag(:) = {Iint};
Mint_real(:) = {Iint};
Mint_imag(:) = {Iint};

%Ordering the elements: PML first
nOfPMLelements = length(inPMLelems);
inINTelems = setdiff(1:nOfElements,inPMLelems);
oElems = [inPMLelems inINTelems];

%Memory allocation for the PML elements
allocation_pml = nOfElementNodes^2*nOfPMLelements;
Ipml = zeros(allocation_pml,1);                              
Jpml = Ipml;                                          
Kxpml_real(:) = {Ipml};
Kxpml_imag(:) = {Ipml};
Kypml_real(:) = {Ipml};
Kypml_imag(:) = {Ipml};
Mpml_real(:)  = {Ipml};
Mpml_imag(:)  = {Ipml};

%Some reshapes, permutes and definitions
N = permute(N,[2 1]);
N = reshape(N,[nOfElementNodes 1 ngauss 1]);        % nOfElementNodes??1??ngauss
Nxi = permute(Nxi,[2 1]);
Nxi = reshape(Nxi,[nOfElementNodes 1 ngauss 1]);    % nOfElementNodes??1??ngauss
Neta = permute(Neta,[2 1]);
Neta = reshape(Neta,[nOfElementNodes 1 ngauss 1]);  % nOfElementNodes??1??ngauss
IPw = reshape(IPw,[1 1 1 1 ngauss]);
np = nOfElementNodes;
nd = 2; %xi,eta
nD = 2; %x,y

%%%% A_I_K = A_i_j_id_jd_g = IPw_g * N_i_id_g * N_j_jd_g
Nxieta = [Nxi Neta];                                % np??nd??ngauss
N_i_id_g = reshape(Nxieta, [np 1 nd 1 ngauss]);
N_i_id_g = bsxfun(@times,IPw,N_i_id_g);
N_j_jd_g = reshape(Nxieta,[1 np 1 nd ngauss]);
A_I_K = bsxfun(@times,N_i_id_g,N_j_jd_g);
A_I_K = reshape(A_I_K,[np*np nd*nd*ngauss]);

%%%% AM_I_K = AM_i_j_g = IPw_g * N_i_g * N_j_g
N_i_g = N;
N_i_g = bsxfun(@times,reshape(IPw,[1 1 ngauss]),N_i_g);
N_j_g = reshape(N,[1,np,ngauss]);
AM_I_K = bsxfun(@times,N_i_g,N_j_g);
AM_I_K = reshape(AM_I_K,[np*np ngauss]);

%Vectorized PML elements
if nOfPMLelements > 0
    vElemsPML = 1:nElemsVectorized:nOfPMLelements;
    if vElemsPML(end) ~= nOfPMLelements
        vElemsPML = [vElemsPML nOfPMLelements];
    end
else 
    vElemsPML = [];
end

%Vectorized INT elements
vElemsINT = nOfPMLelements+1:nElemsVectorized:nOfElements;
if vElemsINT(end) ~= nOfElements
    vElemsINT = [vElemsINT nOfElements];
end

%Total vectorized elements (ordered)
vElems = [vElemsPML(2:end) vElemsINT(2:end)];

%Auxiliar variables
T_t = T';
X_t = X';

%Loop in 2D vectorized elements
iniElem = 1;
for iElem = vElems
    ivElems = oElems(iniElem:iElem);
    nElems = length(ivElems);
    iTe = T_t(:,ivElems);                                                                   % nOfElementNodes??e
    Te = reshape(iTe,1,nElems*nOfElementNodes);                                             % 1??nOfElementNodes*e
    Xe_t = reshape(X_t(:,Te),2,nOfElementNodes,1,nElems);                                   % nD??nOfElementNodes??1??e
    for i = 1:nOfTerms1, coef1e{i} = reshape(vcoef1(Te,i),nOfElementNodes,1,1,nElems); end  % nOfElementNodes??1??1??e
    for i = 1:nOfTerms2, coef2e{i} = reshape(vcoef2(Te,i),nOfElementNodes,1,1,nElems); end  % nOfElementNodes??1??1??e
    for i = 1:nOfTerms3, coef3e{i} = reshape(vcoef3(Te,i),nOfElementNodes,1,1,nElems); end  % nOfElementNodes??1??1??e
    for i = 1:nOfTerms4, coef4e{i} = reshape(vcoef4(Te,i),nOfElementNodes,1,1,nElems); end  % nOfElementNodes??1??1??e
    
    % Assembling parameters
    Te = reshape(Te,1,nOfElementNodes,nElems);
    Te_transp = permute(Te,[2 1 3]);
    aux_row = Te_transp(:,aux_ones,:);
    aux_col = Te(aux_ones,:,:);
    indK = (iniElem-1)*nOfElementNodes^2+1:iElem*nOfElementNodes^2;

    % Elemental Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iElem <= nOfPMLelements %PML elements
        [coef2,coef3,coef4] = computeParameters_PML(N,coef2e,coef3e,coef4e,ngauss,nElems);

        [Mepml,Kexpml,Keypml] = computeElementalMatrices_tensorProduct_PML(...
            Xe_t,coef2,coef3,coef4,Nxieta,np,nd,nD,ngauss,nElems,...
            A_I_K,AM_I_K);
        
        % Assembling for PML area 
        Ipml(indK) = aux_row(:);
        Jpml(indK) = aux_col(:);
        for i = 1:nOfTerms2
            Kxpml_real{i}(indK) = real(Kexpml{i});
            Kxpml_imag{i}(indK) = imag(Kexpml{i});
        end
        for i = 1:nOfTerms3
            Kypml_real{i}(indK) = real(Keypml{i});
            Kypml_imag{i}(indK) = imag(Keypml{i});
        end
        for i = 1:nOfTerms4
            Mpml_real{i}(indK) = real(Mepml{i});
            Mpml_imag{i}(indK) = imag(Mepml{i});
        end
        
    else %Interior elements
        [coef1,coef4] = computeParameters(N,coef1e,coef4e,ngauss,nElems);

        [Keint,Meint] = computeElementalMatrices_tensorProduct(...
            Xe_t,coef1,coef4,Nxieta,np,nd,nD,ngauss,nElems,...
            A_I_K,AM_I_K);
        
        % Assembling for interior area 
        Iint(indK) = aux_row(:);
        Jint(indK) = aux_col(:);
        for i = 1:nOfTerms1
            Kint_real{i}(indK) = real(Keint{i});
            Kint_imag{i}(indK) = imag(Keint{i});
        end
        for i = 1:nOfTerms4
            Mint_real{i}(indK) = real(Meint{i});
            Mint_imag{i}(indK) = imag(Meint{i});
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update
    iniElem = iElem + 1;
end

%Build sparse matrices
pospml = Ipml > 0;
posint = Iint > 0;
for i = 1:nOfTerms1
    Kint{i} = sparse(Iint(posint),Jint(posint),...
        Kint_real{i}(posint) + sqrt(-1)*Kint_imag{i}(posint),...
        nOfNodes,nOfNodes);
    Kint_real{i} = [];
    Kint_imag{i} = []; 
end
for i = 1:nOfTerms2
    Kxpml{i} = sparse(Ipml(pospml),Jpml(pospml),...
        Kxpml_real{i}(pospml) + sqrt(-1)*Kxpml_imag{i}(pospml),...
        nOfNodes,nOfNodes);
    Kxpml_real{i} = [];
    Kxpml_imag{i} = [];
end
for i = 1:nOfTerms3
    Kypml{i} = sparse(Ipml(pospml),Jpml(pospml),...
        Kypml_real{i}(pospml) + sqrt(-1)*Kypml_imag{i}(pospml),...
        nOfNodes,nOfNodes);
    Kypml_real{i} = [];
    Kypml_imag{i} = [];
end
for i = 1:nOfTerms4
    Mint{i} = sparse(Iint(posint),Jint(posint),...
        Mint_real{i}(posint) + sqrt(-1)*Mint_imag{i}(posint),...
        nOfNodes,nOfNodes);
    Mint_real{i} = [];
    Mint_imag{i} = [];
    Mpml{i} = sparse(Ipml(pospml),Jpml(pospml),...
        Mpml_real{i}(pospml) + sqrt(-1)*Mpml_imag{i}(pospml),...
        nOfNodes,nOfNodes);
    Mpml_real{i} = [];
    Mpml_imag{i} = [];
end


function [coef1,coef4] = computeParameters(N,coef1e,coef4e,ngauss,nElems)

auxzeros = zeros(1,1,ngauss,nElems);
nOfTerms1 = numel(coef1e);
nOfTerms4 = numel(coef4e);
coef1 = cell(nOfTerms1,1);
coef4 = cell(nOfTerms4,1);
coef1(:) = {auxzeros};
coef4(:) = {auxzeros};
for i = 1:nOfTerms1, coef1{i} = sum(bsxfun(@times,N,coef1e{i}),1); end   % 1??1??ngauss??e
for i = 1:nOfTerms4, coef4{i} = sum(bsxfun(@times,N,coef4e{i}),1); end   % 1??1??ngauss??e


function [coef2,coef3,coef4] = computeParameters_PML(N,coef2e,coef3e,coef4e,ngauss,nElems)

auxzeros = zeros(1,1,ngauss,nElems);
nOfTerms2 = numel(coef2e);
nOfTerms3 = numel(coef3e);
nOfTerms4 = numel(coef4e);
coef2 = cell(nOfTerms2,1);
coef3 = cell(nOfTerms3,1);
coef4 = cell(nOfTerms4,1);
coef2(:) = {auxzeros};
coef3(:) = {auxzeros};
coef4(:) = {auxzeros};
for i = 1:nOfTerms2, coef2{i} = sum(bsxfun(@times,N,coef2e{i}),1); end   % 1??1??ngauss??e
for i = 1:nOfTerms3, coef3{i} = sum(bsxfun(@times,N,coef3e{i}),1); end   % 1??1??ngauss??e
for i = 1:nOfTerms4, coef4{i} = sum(bsxfun(@times,N,coef4e{i}),1); end   % 1??1??ngauss??e


function [Keint,Meint] = computeElementalMatrices_tensorProduct(...
    Xe_t,coef1,coef4,Nxieta,np,nd,nD,ngauss,nElems,...
    A_I_K,AM_I_K)

%Jacobian
Nxieta_r = reshape(Nxieta,[1 np nd*ngauss]);
Xe_t_p = permute(Xe_t,[1 4 2 3]);
Xe_t_p = reshape(Xe_t_p,[nD*nElems np]);
Nxieta_r_p = permute(Nxieta_r,[2 3 1]);
Nxieta_r_p = reshape(Nxieta_r_p,[np nd*ngauss]);

J = Xe_t_p*Nxieta_r_p;
J = reshape(J,[nD nElems nd ngauss]);
J = permute(J,[1 3 4 2]);

detJ = J(1,1,:,:) .* J(2,2,:,:) - J(1,2,:,:) .* J(2,1,:,:);

invJ11 =  J(2,2,:,:); %Not 1/detJ (optimized)!
invJ12 = -J(1,2,:,:);
invJ21 = -J(2,1,:,:);
invJ22 =  J(1,1,:,:);

invDetJ = 1./detJ;

%Stiffness matrix Keint
nOfTerms = numel(coef1);
Keint = cell(nOfTerms,1);
Keint(:) = {zeros(np*np,nElems)};
B_K = zeros(nd,nd,ngauss,nElems);
B_K(1,1,:,:) = invDetJ.*(invJ11.*invJ11 + invJ12.*invJ12);
B_K(1,2,:,:) = invDetJ.*(invJ11.*invJ21 + invJ12.*invJ22);
B_K(2,1,:,:) = invDetJ.*(invJ21.*invJ11 + invJ22.*invJ12);
B_K(2,2,:,:) = invDetJ.*(invJ21.*invJ21 + invJ22.*invJ22);
for i = 1:nOfTerms
    iB_K = bsxfun(@times,B_K,coef1{i});

    iB_K = reshape(iB_K,[nd*nd*ngauss nElems]);
    Keint{i} = A_I_K*iB_K;
    Keint{i} = reshape(Keint{i},[np np nElems]);
end

%Mass matrix Mpml
nOfTerms = numel(coef4);
Meint = cell(nOfTerms,1);
Meint(:) = {zeros(np*np,nElems)};
for i = 1:nOfTerms
    Meint{i} = AM_I_K*reshape(detJ.*coef4{i},[ngauss nElems]);
    Meint{i} = reshape(Meint{i},[np np nElems]);
end


function [Mepml,Kexpml,Keypml] = computeElementalMatrices_tensorProduct_PML(...
    Xe_t,coef2,coef3,coef4,Nxieta,np,nd,nD,ngauss,nElems,...
    A_I_K,AM_I_K)

%Jacobian
Nxieta_r = reshape(Nxieta,[1 np nd*ngauss]);
Xe_t_p = permute(Xe_t,[1 4 2 3]);
Xe_t_p = reshape(Xe_t_p,[nD*nElems np]);
Nxieta_r_p = permute(Nxieta_r,[2 3 1]);
Nxieta_r_p = reshape(Nxieta_r_p,[np nd*ngauss]);

J = Xe_t_p*Nxieta_r_p;
J = reshape(J,[nD nElems nd ngauss]);
J = permute(J,[1 3 4 2]);

detJ = J(1,1,:,:) .* J(2,2,:,:) - J(1,2,:,:) .* J(2,1,:,:);

invJ11 =  J(2,2,:,:); %Not 1/detJ (optimized)!
invJ12 = -J(1,2,:,:);
invJ21 = -J(2,1,:,:);
invJ22 =  J(1,1,:,:);

invDetJ = 1./detJ;

invJ_11_11 = invDetJ.*invJ11.*invJ11;
invJ_12_12 = invDetJ.*invJ12.*invJ12;
invJ_11_21 = invDetJ.*invJ11.*invJ21;
invJ_12_22 = invDetJ.*invJ12.*invJ22;
invJ_21_11 = invDetJ.*invJ21.*invJ11;
invJ_22_12 = invDetJ.*invJ22.*invJ12;
invJ_21_21 = invDetJ.*invJ21.*invJ21;
invJ_22_22 = invDetJ.*invJ22.*invJ22;

%Stiffness matrix Kxpml
nOfTerms = numel(coef2);
Kexpml = cell(nOfTerms,1);
Kexpml(:) = {zeros(np*np,nElems)};
B_K = zeros(nd,nd,ngauss,nElems);
B_K(1,1,:,:) = invJ_11_11;
B_K(1,2,:,:) = invJ_11_21;
B_K(2,1,:,:) = invJ_21_11;
B_K(2,2,:,:) = invJ_21_21;
for i = 1:nOfTerms
    iB_K = bsxfun(@times,B_K,coef2{i});

    iB_K = reshape(iB_K,[nd*nd*ngauss nElems]);
    Kexpml{i} = A_I_K*iB_K;
    Kexpml{i} = reshape(Kexpml{i},[np np nElems]);
end

%Stiffness matrix Kypml
nOfTerms = numel(coef3);
Keypml = cell(nOfTerms,1);
Keypml(:) = {zeros(np*np,nElems)};
B_K = zeros(nd,nd,ngauss,nElems);
B_K(1,1,:,:) = invJ_12_12;
B_K(1,2,:,:) = invJ_12_22;
B_K(2,1,:,:) = invJ_22_12;
B_K(2,2,:,:) = invJ_22_22;
for i = 1:nOfTerms
    iB_K = bsxfun(@times,B_K,coef3{i});

    iB_K = reshape(iB_K,[nd*nd*ngauss nElems]);
    Keypml{i} = A_I_K*iB_K;
    Keypml{i} = reshape(Keypml{i},[np np nElems]);
end

%Mass matrix Mpml
nOfTerms = numel(coef4);
Mepml = cell(nOfTerms,1);
Mepml(:) = {zeros(np*np,nElems)};
for i = 1:nOfTerms
    Mepml{i} = AM_I_K*reshape(detJ.*coef4{i},[ngauss nElems]);
    Mepml{i} = reshape(Mepml{i},[np np nElems]);
end

















