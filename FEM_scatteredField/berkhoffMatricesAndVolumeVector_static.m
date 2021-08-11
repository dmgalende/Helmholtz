function [K,f,elementalInfo,newPositionFaceNodes] = berkhoffMatricesAndVolumeVector_static...
    (X,T,theReferenceElement,amp0,kvector,kVec,ccgVec,sigmaVec,omega,NNZ)

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
%  NNZ: estimated number of non zero entries of the matrix (NON NECESSARY
%       with this version of the function)
%  Output:
%  K: (mass - stiffness) matrix ---ONLY FACE NODES---
%  f: volume vector ---ONLY FACE NODES---
%  elementalInfo: info needed for postprocess the inner nodes solution
%  newPositionFaceNodes: new global numbering for each face node which is
%  going to be used to solve the system (zero numbering is used for inner
%  nodes)

%Reference element numbering
innerNodes = theReferenceElement.innerNodes;
nOfInnerNodes = length(innerNodes);
faceNodes = unique(theReferenceElement.faceNodes);
nOfFaceNodes = length(faceNodes);

%Number of elements and number of mesh nodes
nOfElements = size(T,1);
nOfTotalInnerNodes = nOfElements*nOfInnerNodes;
nOfTotalNodes = size(X,1);
nOfNodes = nOfTotalNodes - nOfTotalInnerNodes;

%Information of the reference element
IPw = theReferenceElement.IPweights;                % ng·1
N = theReferenceElement.N;                          % ng·nOfElementNodes
Nxi = theReferenceElement.Nxi;                      % ng·nOfElementNodes
Neta = theReferenceElement.Neta;                    % ng·nOfElementNodes
ngauss = length(IPw);

%Inner nodes position in global matrix (sorted)
innerNodesBool = false(1,nOfTotalNodes);
innerNodesBool(sort(reshape(T(:,innerNodes),nOfTotalInnerNodes,1),'ascend')) = true;

%New connectivity for face nodes (same sparsivity)
newPositionFaceNodes = cumsum(~innerNodesBool);
newPositionFaceNodes(innerNodesBool) = 0;

%Memory allocation
nOfElementNodes = size(T,2);
nOfElementNodesMod = nOfElementNodes - nOfInnerNodes;
allocationF = nOfElementNodesMod*nOfElements;
f = zeros(allocationF,1);
IF = f;
aux_ones = ones(1,nOfElementNodesMod);
allocation = nOfElementNodesMod^2*nOfElements;
I = zeros(allocation,1);                            % allocation·1
J = I;                                              % allocation·1
Kreal = I;                                          % allocation·1
Kimag = I;                                          % allocation·1

%Some reshapes and permutes
N = permute(N,[2 1]);
N = reshape(N,[nOfElementNodes 1 ngauss 1]);        % nOfElementNodes·1·ngauss
Nxi = permute(Nxi,[2 1]);
Nxi = reshape(Nxi,[nOfElementNodes 1 ngauss 1]);    % nOfElementNodes·1·ngauss
Neta = permute(Neta,[2 1]);
Neta = reshape(Neta,[nOfElementNodes 1 ngauss 1]);  % nOfElementNodes·1·ngauss
IPw = reshape(IPw,[1 1 1 1 ngauss]);

%%%% A_I_K = A_i_j_id_jd_g = IPw_g * N_i_id_g * N_j_jd_g
np = nOfElementNodes;
nd = 2; %xi,eta
nD = 2; %x,y

Nxieta = [Nxi Neta];                                % np·nd·ngauss
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

%%%% AF1_I_K = AF1_i_j_jd_g = IPw_g * N_i_g * N_j_jd_g
N_j_jd_g = reshape(Nxieta,[1 np nd ngauss]);
AF1_I_K = bsxfun(@times,reshape(N_i_g,[np 1 1 ngauss]),N_j_jd_g);
AF1_I_K = reshape(AF1_I_K,[np nd*np*ngauss]);

%%%% AF2_I_K = AF2_i_g = IPw_g * N_i_g
AF2_I_K = reshape(N_i_g,[np ngauss]);

%Vectorized elements
nElemsVectorized = 1024;
vElems = 1:nElemsVectorized:nOfElements;
if vElems(end) ~= nOfElements
    vElems = [vElems nOfElements];
end
vElems(1) = [];
iniElem = 1;

%Auxiliar variables
T_t = T';
X_t = X';
sigmaVec_t = sigmaVec.';

%Elemental info
elementalInfo = zeros(nOfInnerNodes,nOfFaceNodes+1,nOfElements);

%Loop in 2D vectorized elements
for iElem = vElems
    ivElems = iniElem:iElem;
    nElems = length(ivElems);
    iTe = T_t(:,ivElems);                                           % nOfElementNodes·e
    Te = reshape(iTe,1,nElems*nOfElementNodes);                     % 1·nOfElementNodes*e
    Xe_t = reshape(X_t(:,Te),2,nOfElementNodes,1,nElems);           % nD·nOfElementNodes·1·e
    Xe = permute(Xe_t,[2 1 3 4]);                                   % nOfElementNodes·nD·1·e
    ccge = reshape(ccgVec(Te),nOfElementNodes,1,1,nElems);          % nOfElementNodes·1·1·e
    ke = reshape(kVec(Te),nOfElementNodes,1,1,nElems);              % nOfElementNodes·1·1·e
    sigmae = reshape(sigmaVec_t(:,Te),2,nOfElementNodes,1,nElems);  % nD·nOfElementNodes·1·e
    sigmae = permute(sigmae,[2 1 3 4]);                             % nOfElementNodes·nD·1·e

    % Elemental Matrix and vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [k,ccg,phi0,phi0der,phi0Lap,param] = computeParameters(...
        Xe,N,ccge,ke,sigmae,amp0,kvector,omega);

    [Ke,fe] = computeMeMinusKe_tensorProduct(...
        Xe_t,k,ccge,param,ccg,phi0,phi0der,phi0Lap,...
        Nxieta,np,nd,nD,ngauss,nElems,A_I_K,AM_I_K,AF1_I_K,AF2_I_K);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Static condensation
    Kes = zeros(nOfFaceNodes,nOfFaceNodes,nElems);
    fes = zeros(nOfFaceNodes,nElems);
    for jElem = 1:nElems
        
        %Variables
        As = Ke(faceNodes,faceNodes,jElem);
        invBs = inv(Ke(innerNodes,innerNodes,jElem)); %inverse!
        cs = Ke(faceNodes,innerNodes,jElem);
        fs = fe(faceNodes,jElem);
        gs = fe(innerNodes,jElem);
        
        %New elemental matrix and vector (only face nodes)
        invBs_csT = invBs*cs.';
        invBs_gs = invBs*gs;
        Kes(:,:,jElem) = As - cs*invBs_csT;
        fes(:,jElem) = fs - cs*invBs_gs;
        
        %Store elemental info for postprocess the solution at inner nodes
        elementalInfo(:,:,ivElems(jElem)) = [invBs_gs invBs_csT];
    end
    
    %New connectivity for face nodes due to static condensation
    Te = newPositionFaceNodes(reshape(iTe(faceNodes,:),1,nElems*nOfElementNodesMod));

    % Assembling
    indF = (iniElem-1)*nOfElementNodesMod+1:iElem*nOfElementNodesMod;
    IF(indF) = Te;
    f(indF) = f(indF) + fes(:);
    Te = reshape(Te,1,nOfElementNodesMod,nElems);
    Te_transp = permute(Te,[2 1 3]);
    aux_row = Te_transp(:,aux_ones,:);
    aux_col = Te(aux_ones,:,:);
    indK = (iniElem-1)*nOfElementNodesMod^2+1:iElem*nOfElementNodesMod^2;
    I(indK) = aux_row(:);
    J(indK) = aux_col(:);
    Kreal(indK) = real(Kes);
    Kimag(indK) = imag(Kes);

    % Update
    iniElem = iElem + 1;
end

f = accumarray(IF,f);
K = sparse(I,J,Kreal + sqrt(-1)*Kimag,nOfNodes,nOfNodes);


function [k,ccg,phi0,phi0der,phi0Lap,param] = computeParameters(...
    Xe,N,ccge,ke,sigmae,amp0,kvector,omega)

ccg = sum(bsxfun(@times,N,ccge),1);                         % 1·1·ngauss·e
k = sum(bsxfun(@times,N,ke),1);                             % 1·1·ngauss·e

xy = sum(bsxfun(@times,N,Xe),1);                            % 1·2·ngauss·e
phi0 = amp0*exp(sqrt(-1)*sum(bsxfun(@times,kvector,xy),2)); % 1·1·ngauss·e
phi0der = sqrt(-1)*bsxfun(@times,kvector,phi0);             % 1·2·ngauss·e
phi0Lap = -sum(kvector.*kvector,2)*phi0;                    % 1·1·ngauss·e

sigma = sum(bsxfun(@times,N,sigmae),1);                     % 1·2·ngauss·e
param = 1 + sqrt(-1)*sigma/omega;                           % 1·2·ngauss·e


function [C,f] = computeMeMinusKe_tensorProduct(...
    Xe_t,k,ccge,param,ccg,phi0,phi0der,phi0Lap,...
    Nxieta,np,nd,nD,ngauss,nElems,A_I_K,AM_I_K,AF1_I_K,AF2_I_K)

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
beta = zeros(nd,ngauss,nElems);
beta(1,:,:) = param1param2.*(invJ11.*phi0der(:,1,:,:) + invJ12.*phi0der(:,2,:,:));
beta(2,:,:) = param1param2.*(invJ21.*phi0der(:,1,:,:) + invJ22.*phi0der(:,2,:,:));
BF1_K = bsxfun(@times,reshape(beta,[1 nd ngauss nElems]),reshape(ccge,[np 1 1 nElems]));
BF1_K = reshape(BF1_K,[np*nd*ngauss nElems]);

BF2_K = detJ.*param1param2.*(ccg.*phi0Lap + k2ccg.*phi0);
BF2_K = reshape(BF2_K,[ngauss nElems]);

fe1 = AF1_I_K*BF1_K;
fe2 = AF2_I_K*BF2_K;
f = fe1 + fe2;

























