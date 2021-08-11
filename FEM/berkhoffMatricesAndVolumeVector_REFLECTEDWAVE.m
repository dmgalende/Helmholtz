function [K,ftotal] = berkhoffMatricesAndVolumeVector...
    (X,T,theReferenceElement,u0Vec,kVec,ccgVec,sigmaVec,omega,inPMLelems,nOfBoundaryElements)

% Input:
%  X: nodal coordinates
%  T: connectivity matrix for 2D elements
%  theReferenceElement: struct with the information of the reference
%  element
%  amp0: amplitude of incident potential
%  kvector: k0*d0, where k0, d0 are the wave number and the direction of
%  incident potential
%  kVec: nodal values of wave number
%  ccgVec: nodal values of the c*cg parameter
%  sigmaVec: PML parameter
%  omega: angular frequency
% Output:
%  K: (mass - stiffness) matrix
%  f: volume vector

%Vectorized elements
nElemsVectorized = 1024;

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfNodes = size(X,1);

%Information of the reference element
IPw = theReferenceElement.IPweights;                % ng??1
N = theReferenceElement.N;                          % ng??nOfElementNodes
Nxi = theReferenceElement.Nxi;                      % ng??nOfElementNodes
Neta = theReferenceElement.Neta;                    % ng??nOfElementNodes
ngauss = length(IPw);

%Memory allocation
nOfElementNodes = size(T,2);
if nargin == 11 %NEFEM computation
    nOfElements_F = nOfElements + nOfBoundaryElements;
else
    nOfElements_F = nOfElements;
end
allocationF = nOfElementNodes*nOfElements_F;
f = zeros(allocationF,1);
ftotal = zeros(nOfNodes,1);
IF = f;
aux_ones = ones(1,nOfElementNodes);
allocation = nOfElementNodes^2*nOfElements;
I = zeros(allocation,1);                            % allocation??1
J = I;                                              % allocation??1
Kreal = I;                                          % allocation??1
Kimag = I;                                          % allocation??1

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

%%%% AF1_I_K = AF1_i_j_id_jd_g = IPw_g * N_i_id_g * N_j_jd_g
AF1_I_K = reshape(A_I_K,[np np*nd*nd*ngauss]);

%%%% AF2_I_K = AF2_i_g = IPw_g * N_i_g
AF2_I_K = reshape(N_i_g,[np ngauss]);

%Ordering the elements: PML first
nOfPMLelements = length(inPMLelems);
inINTelems = setdiff(1:nOfElements,inPMLelems);
oElems = [inPMLelems inINTelems];

%Vectorized PML elements
vElemsPML = 1:nElemsVectorized:nOfPMLelements;
if vElemsPML(end) ~= nOfPMLelements
    vElemsPML = [vElemsPML nOfPMLelements];
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
sigmaVec_t = sigmaVec.';

%Loop in 2D vectorized elements
iniElem = 1;
for iElem = vElems
    ivElems = oElems(iniElem:iElem);
    nElems = length(ivElems);
    iTe = T_t(:,ivElems);                                           % nOfElementNodes??e
    Te = reshape(iTe,1,nElems*nOfElementNodes);                     % 1??nOfElementNodes*e
    Xe_t = reshape(X_t(:,Te),2,nOfElementNodes,1,nElems);           % nD??nOfElementNodes??1??e
    ccge = reshape(ccgVec(Te),nOfElementNodes,1,1,nElems);          % nOfElementNodes??1??1??e
    ke = reshape(kVec(Te),nOfElementNodes,1,1,nElems);              % nOfElementNodes??1??1??e
    u0e = reshape(u0Vec(Te),nOfElementNodes,1,1,nElems);            % nOfElementNodes??1??1??e
    sigmae = reshape(sigmaVec_t(:,Te),2,nOfElementNodes,1,nElems);  % nD??nOfElementNodes??1??e
    sigmae = permute(sigmae,[2 1 3 4]);                             % nOfElementNodes??nD??1??e

    % Elemental Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iElem <= nOfPMLelements %PML elements --> where f=0
        [k,ccg,param] = computeParameters_PML(N,ccge,ke,sigmae,omega);

        [Ke,fe] = computeMeMinusKe_tensorProduct_PML(...
            Xe_t,k,param,ccg,Nxieta,np,nd,nD,ngauss,nElems,...
            A_I_K,AM_I_K);
    else %Interior elements --> where f~=0 and param=1
        [k,ccg,u0] = computeParameters(N,ccge,ke,u0e);

        [Ke,fe] = computeMeMinusKe_tensorProduct(...
            Xe_t,k,ccg,u0e,u0,Nxieta,np,nd,nD,ngauss,nElems,...
            A_I_K,AM_I_K,AF1_I_K,AF2_I_K);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Assembling
    indF = (iniElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    IF(indF) = Te;
    f(indF) = f(indF) + fe(:);
    Te = reshape(Te,1,nOfElementNodes,nElems);
    Te_transp = permute(Te,[2 1 3]);
    aux_row = Te_transp(:,aux_ones,:);
    aux_col = Te(aux_ones,:,:);
    indK = (iniElem-1)*nOfElementNodes^2+1:iElem*nOfElementNodes^2;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Kreal(indK) = real(Ke);
    Kimag(indK) = imag(Ke);

    % Update
    iniElem = iElem + 1;
end

%Vector and sparse matrix
pos0 = find(IF == 0);
IF(pos0) = pos0; %Trick to add zero values. This match dimensions in nefem computations
f = accumarray(IF,f);
nodesInterior = unique(T);
ftotal(nodesInterior) = f(nodesInterior);
K = sparse(I,J,Kreal + sqrt(-1)*Kimag,nOfNodes,nOfNodes);

%Correction for vector ftotal in PML
% nodesPML = unique(T(inPMLelems,:));
% ftotal(nodesPML) = 0;


function [k,ccg,u0] = computeParameters(N,ccge,ke,u0e)

ccg = sum(bsxfun(@times,N,ccge),1);                         % 1??1??ngauss??e
k = sum(bsxfun(@times,N,ke),1);                             % 1??1??ngauss??e
u0 = sum(bsxfun(@times,N,u0e),1);                           % 1??1??ngauss??e

function [k,ccg,param] = computeParameters_PML(N,ccge,ke,sigmae,omega)

ccg = sum(bsxfun(@times,N,ccge),1);                         % 1??1??ngauss??e
k = sum(bsxfun(@times,N,ke),1);                             % 1??1??ngauss??e
sigma = sum(bsxfun(@times,N,sigmae),1);                     % 1??2??ngauss??e
param = 1 + sqrt(-1)*sigma/omega;                           % 1??2??ngauss??e


function [C,f] = computeMeMinusKe_tensorProduct(...
    Xe_t,k,ccg,u0e,u0,Nxieta,np,nd,nD,ngauss,nElems,...
    A_I_K,AM_I_K,AF1_I_K,AF2_I_K)

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

%Tensor for stiffness matrix
B_K = zeros(nd,nd,ngauss,nElems);
B_K(1,1,:,:) = invDetJ.*(invJ11.*invJ11 + invJ12.*invJ12);
B_K(1,2,:,:) = invDetJ.*(invJ11.*invJ21 + invJ12.*invJ22);
B_K(2,1,:,:) = invDetJ.*(invJ21.*invJ11 + invJ22.*invJ12);
B_K(2,2,:,:) = invDetJ.*(invJ21.*invJ21 + invJ22.*invJ22);
B_K = bsxfun(@times,B_K,ccg);

B_K = reshape(B_K,[nd*nd*ngauss nElems]);
Ke = A_I_K*B_K;
Ke = reshape(Ke,[np np nElems]);

%Tensor for mass matrix
detJk2ccg = detJ.*k.^2.*ccg;
BM_K = reshape(detJk2ccg,[ngauss nElems]);

Me = AM_I_K*BM_K;
Me = reshape(Me,[np np nElems]);

%Mass minus stiffness matrix
C = Me - Ke;

%Independent vector
BF1_K = bsxfun(@times,reshape(B_K,[1 nd*nd*ngauss nElems]),reshape(u0e,[np 1 nElems]));
BF1_K = reshape(BF1_K,[np*nd*nd*ngauss nElems]);

BF2_K = detJk2ccg.*u0;
BF2_K = reshape(BF2_K,[ngauss nElems]);

fe1 = AF1_I_K*BF1_K;
fe2 = AF2_I_K*BF2_K;
f = fe2 - fe1;

function [C,f] = computeMeMinusKe_tensorProduct_PML(...
    Xe_t,k,param,ccg,Nxieta,np,nd,nD,ngauss,nElems,...
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

%Tensor for stiffness matrix
alpha = zeros(1,nD,ngauss,nElems);
alpha(:,1,:,:) = ccg.*param(1,2,:,:)./param(1,1,:,:);
alpha(:,2,:,:) = ccg.*param(1,1,:,:)./param(1,2,:,:);

B_K = zeros(nd,nd,ngauss,nElems);
B_K(1,1,:,:) = invDetJ.*(alpha(1,1,:,:).*invJ11.*invJ11 + alpha(1,2,:,:).*invJ12.*invJ12);
B_K(1,2,:,:) = invDetJ.*(alpha(1,1,:,:).*invJ11.*invJ21 + alpha(1,2,:,:).*invJ12.*invJ22);
B_K(2,1,:,:) = invDetJ.*(alpha(1,1,:,:).*invJ21.*invJ11 + alpha(1,2,:,:).*invJ22.*invJ12);
B_K(2,2,:,:) = invDetJ.*(alpha(1,1,:,:).*invJ21.*invJ21 + alpha(1,2,:,:).*invJ22.*invJ22);

B_K = reshape(B_K,[nd*nd*ngauss nElems]);
Ke = A_I_K*B_K;
Ke = reshape(Ke,[np np nElems]);

%Tensor for mass matrix
param1param2 = param(:,1,:,:).*param(:,2,:,:);
k2ccg = k.^2.*ccg;
BM_K = detJ.*k2ccg.*param1param2;
BM_K = reshape(BM_K,[ngauss nElems]);

Me = AM_I_K*BM_K;
Me = reshape(Me,[np np nElems]);

%Mass minus stiffness matrix
C = Me - Ke;

%Independent vector
f = zeros(np,nElems);

















