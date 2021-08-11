function [KminusM,f,R] = berkhoffMatricesAndVolumeVector_CDG_P_adapt...
    (X,T,theReferenceElement,A,kvector,k,alpha,sigma,omega,p_adaptivity)

nOfElements = size(T,1);
elements = 1:nOfElements;
elem_p = p_adaptivity.elem_p;
elem_check = p_adaptivity.elem_check;

% Allocation
allocation = 0.25*sum(((elem_p+1).*(elem_p+2)).^2);
dim = 0.5*sum((elem_p+1).*(elem_p+2));
I = zeros(allocation,1);
J = zeros(allocation,1);
Kr = zeros(allocation,1);
Ki = zeros(allocation,1);
f = zeros(0.5*sum((elem_p+1).*(elem_p+2)),1);
R = cell(nOfElements,1);

% number of nodes in the base elements
nOfElementNodes = size(theReferenceElement.NodesCoord,1);

% loop in base elements (not refined)
COUNT = 0.1;
time = cputime;
loop_count = 0;
nOfBaseElements = numel(elements(~elem_check));
disp(' ')
disp(['Loop in ' num2str(nOfBaseElements) ' base elements'])
for iElem = elements(~elem_check)
    loop_count = loop_count+1;
    Te = T(iElem,:);
    alpha_el = alpha(Te);
    ke = k(Te);
    sigma_el = sigma(Te,:);
    Xe = X(Te,:);

    % elemental matrices
    [KeminusMe,fe,Re] = KeMeElementalMatrices(Xe,theReferenceElement,A,ke,kvector,alpha_el,sigma_el,omega,nOfElementNodes);

    % NEW ASSEMBLY
    ind = 0.5*sum((elem_p(1:iElem-1)+1).*(elem_p(1:iElem-1)+2))+1 : ...
        0.5*sum((elem_p(1:iElem)+1).*(elem_p(1:iElem)+2));
    ind_transp = transpose(ind);
    aux_row = ind_transp(:,ones(1,nOfElementNodes));
    aux_col = ind(ones(1,nOfElementNodes),:);
    indK = 0.25*sum(((elem_p(1:iElem-1)+1).*(elem_p(1:iElem-1)+2)).^2)+1 : ...
        0.25*sum(((elem_p(1:iElem)+1).*(elem_p(1:iElem)+2)).^2);
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Kr(indK) = real(KeminusMe);
    Ki(indK) = imag(KeminusMe);
    f(ind) = f(ind) + fe;
    R{iElem} = sparse(Re);
    clear KeminusMe fe Re
    if loop_count > nOfBaseElements*COUNT
        disp(['Calculated: ' num2str(100*COUNT) '% ; time elapsed ' num2str(cputime-time)])
        time = cputime;
        COUNT = COUNT + 0.1;
    end
end

% loop in refined elements
COUNT = 0.1;
time = cputime;
loop_count = 0;
nOfRefinedElements = numel(elements(elem_check));
disp(' ')
disp(['Loop in ' num2str(nOfRefinedElements) ' refined elements'])
for iElem = elements(elem_check)
    loop_count = loop_count+1;
    p = elem_p(iElem);
    refEl = p_adaptivity.referenceElements.(['P' num2str(p)]);
    nOfElementNodes = size(refEl.NodesCoord,1);

    % parameters on the nodes of the refined elements
    indRef = 0.5*sum((elem_p(elem_check(1:iElem-1))+1).*...
        (elem_p(elem_check(1:iElem-1))+2))+(1:nOfElementNodes);
    Xe = p_adaptivity.coord(indRef,:);
    alpha_el = p_adaptivity.ccg(indRef);
    ke = p_adaptivity.waveNumber(indRef);
    sigma_el = p_adaptivity.sigma(indRef,:);

    % elemental matrices
    [KeminusMe,fe,Re] = KeMeElementalMatrices(Xe,refEl,A,ke,kvector,alpha_el,sigma_el,omega,nOfElementNodes);

    % NEW ASSEMBLY
    ind = 0.5*sum((elem_p(1:iElem-1)+1).*(elem_p(1:iElem-1)+2))+1 : ...
        0.5*sum((elem_p(1:iElem)+1).*(elem_p(1:iElem)+2));
    ind_transp = transpose(ind);
    aux_row = ind_transp(:,ones(1,nOfElementNodes));
    aux_col = ind(ones(1,nOfElementNodes),:);
    indK = 0.25*sum(((elem_p(1:iElem-1)+1).*(elem_p(1:iElem-1)+2)).^2)+1 : ...
        0.25*sum(((elem_p(1:iElem)+1).*(elem_p(1:iElem)+2)).^2);
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Kr(indK) = real(KeminusMe);
    Ki(indK) = imag(KeminusMe);
    f(ind) = f(ind) + fe;
    R{iElem} = sparse(Re);
    clear KeminusMe fe Re
    if loop_count > nOfRefinedElements*COUNT
        disp(['Calculated: ' num2str(100*COUNT) '% ; time elapsed ' num2str(cputime-time)])
        time = cputime;
        COUNT = COUNT + 0.1;
    end
end
disp(' ')
disp('Creating sparse matrix')
disp(' ')
KminusM = sparse(I(I~=0),J(I~=0),Kr(I~=0)+sqrt(-1)*Ki(I~=0),dim,dim);

function [KeminusMe,fe,Re] = KeMeElementalMatrices(Xe,theReferenceElement,A,ke,kvector,...
    ccge,sigma,omega,nOfElementNodes)

KeminusMe = zeros(nOfElementNodes,nOfElementNodes);
fe = zeros(nOfElementNodes,1);
Me_CCg_P = zeros(2*nOfElementNodes,2*nOfElementNodes);
Me = Me_CCg_P;

%Information of the reference element
IPw = theReferenceElement.IPweights;
N = theReferenceElement.N;
Nxi = theReferenceElement.Nxi;
Neta = theReferenceElement.Neta;

% NN = expandMatrix(N);
% the following lines replace expandMatrix
nn =2;
NN = zeros([size(N) nn nn]);
NN(:,:,1:nn+1:nn^2) = repmat(N, [1 1 nn]);
NN = permute(NN, [3 1 4 2]);
NN = reshape(NN, nn*size(N));

%Number of Gauss points
ngauss = length(IPw);

% x and y coordinates of the element nodes
xe = Xe(:,1); ye = Xe(:,2);

%LOOP IN GAUSS POINTS
for g = 1:ngauss

    %Shape functions and derivatives at the current integration point
    N_g = N(g,:);
    NN_g = NN([2*g-1 2*g],:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);

    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye
        Neta_g*xe  Neta_g*ye];
    if det(J)<0
        error('computeElementalMatrices: det(J)<0')
    end

    %Integration weight
    dvolu=IPw(g)*det(J);

    %x and y derivatives
    invJ = inv(J);
    Nx_g = invJ(1,1)*Nxi_g + invJ(1,2)*Neta_g;
    Ny_g = invJ(2,1)*Nxi_g + invJ(2,2)*Neta_g;

    %phase and group celerities at current gauss point
    ccg = N_g*ccge;
    ccgx = Nx_g*ccge;
    ccgy = Ny_g*ccge;
    ccgGrad = [ccgx ccgy];
    k = N_g*ke;

    %PML parameters
    sigma_g = N_g*sigma;
    param = 1 + sqrt(-1)*sigma_g/omega;

    %Incident potential and derivative at the current integration point
    xy_g = N_g*Xe;
    phi0_g = A*exp(i*kvector*xy_g');
    phi0der_g = i*kvector'*phi0_g;
    phi0Lap_g = -kvector*kvector'*phi0_g;
    P_NN_g(1,:) = NN_g(1,:)*param(2)/param(1);
    P_NN_g(2,:) = NN_g(2,:)*param(1)/param(2);

    %Contribution of the current integration point to the elemental matrix
    KeminusMe = KeminusMe + ccg*(param(2)/param(1)*Nx_g'*Nx_g +...
        param(1)/param(2)*Ny_g'*Ny_g - k^2*param(1)*param(2)*N_g'*N_g)*dvolu;
    fe = fe + N_g'*param(1)*param(2)*(ccgGrad*phi0der_g + ccg*phi0Lap_g + ...
        k^2*ccg*phi0_g)*dvolu;
    Me_CCg_P = Me_CCg_P + ccg*transpose(P_NN_g)*NN_g*dvolu;
    Me = Me + NN_g'*NN_g*dvolu;
end
inv_Me = inv(Me);
Re = inv_Me*Me_CCg_P*inv_Me;
