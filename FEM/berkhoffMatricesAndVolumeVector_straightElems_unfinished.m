function [K,ftotal] = berkhoffMatricesAndVolumeVector_straightElems...
    (X,T,theReferenceElement,amp0,kvector,kVec,ccgVec,sigmaVec,omega,nOfBoundaryElements)

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

%Some reshapes and permutes
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

%%%% AMS_I_K = AMS_k_i_j = sum_g{IPw_g * N_k_g * N_i_g * N_j_g}
N_i_g = N;
N_i_g = bsxfun(@times,reshape(IPw,[1 1 ngauss]),N_i_g);
AMS_I_K = bsxfun(@times,reshape(N_i_g,[np 1 1 ngauss]),...
                 reshape(N,[1 np 1 ngauss]));
AMS_I_K = bsxfun(@times,AMS_I_K,reshape(N,[1 1 np ngauss]));
AMS_I_K = sum(AMS_I_K,4);
AMS_I_K = reshape(AMS_I_K,[np*np np]);

%%%% AKS1_I_K = AKS1_k_i_j = sum_g{IPw_g * Nxi_i_g * Nxi_j_g * N_k_g}
AKS_xi_xi_aux = bsxfun(@times,reshape(N_i_g,[np 1 1 ngauss]),...
                 reshape(Nxi,[1 np 1 ngauss]));
AKS_xi_xi = bsxfun(@times,AKS_xi_xi_aux,reshape(Nxi,[1 1 np ngauss]));
AKS_xi_xi = sum(AKS_xi_xi,4);
AKS_xi_xi = permute(AKS_xi_xi,[2 3 1]);
AKS_xi_xi = reshape(AKS_xi_xi,[np*np np]);

%%%% AKS2_I_K = AKS2_k_i_j = sum_g{IPw_g * Nxi_i_g * Neta_j_g * N_k_g}
AKS_xi_eta_aux = bsxfun(@times,AKS_xi_xi_aux,reshape(Neta,[1 1 np ngauss]));
AKS_xi_eta_aux = sum(AKS_xi_eta_aux,4);
AKS_xi_eta_aux = permute(AKS_xi_eta_aux,[2 3 1]);
AKS_xi_eta = reshape(AKS_xi_eta_aux,[np*np np]);

%%%% AKS3_I_K = AKS3_k_i_j = sum_g{IPw_g * Neta_i_g * Nxi_j_g * N_k_g}
AKS_eta_xi = permute(AKS_xi_eta_aux,[2 1 3]);
AKS_eta_xi = reshape(AKS_eta_xi,[np*np np]);
AKS_sym = AKS_xi_eta + AKS_eta_xi;

%%%% AKS4_I_K = AKS4_k_i_j = sum_g{IPw_g * N_k_g * Neta_i_g * Neta_j_g}
AKS_eta_eta = bsxfun(@times,reshape(N_i_g,[np 1 1 ngauss]),...
                 reshape(Neta,[1 np 1 ngauss]));
AKS_eta_eta = bsxfun(@times,AKS_eta_eta,reshape(Neta,[1 1 np ngauss]));
AKS_eta_eta = sum(AKS_eta_eta,4);
AKS_eta_eta = permute(AKS_eta_eta,[2 3 1]);
AKS_eta_eta = reshape(AKS_eta_eta,[np*np np]);

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

%Loop in 2D vectorized elements
for iElem = vElems
    ivElems = iniElem:iElem;
    nElems = length(ivElems);
    iTe = T_t(:,ivElems);                                           % nOfElementNodes??e
    Te = reshape(iTe,1,nElems*nOfElementNodes);                     % 1??nOfElementNodes*e
    
    %Parameters
    ke = reshape(kVec(Te),np,nElems);                               % nOfElementNodes??e
    ccge = reshape(ccgVec(Te),np,nElems);
    sigmae = reshape(sigmaVec_t(:,Te),2,np,nElems);
    parame = 1 + sqrt(-1)*sigmae/omega;
    
    %Constant jacobian
    Xe_t = reshape(X_t(:,Te),2,np,nElems);                          % nD??nOfElementNodes??e
    Xe = permute(Xe_t,[2 1 3]);                                     % nOfElementNodes??nD??e
    v1 = Xe(1,:,:); v2 = Xe(2,:,:); v3 = Xe(3,:,:);
    J1 = reshape((v2 - v1)/2,nD,nElems); 
    J2 = reshape((v3 - v1)/2,nD,nElems);
    detJ = J1(1,:) .* J2(2,:) - J2(1,:) .* J1(2,:);                 % 1??e
    invJ11 = J2(2,:); invJ12 = -J1(2,:);
    invJ21 = -J2(1,:); invJ22 = J1(1,:);
    invdetJ = 1./detJ;
    
    %Mass 
    param1 = reshape(parame(1,:,:),[np nElems]);
    param2 = reshape(parame(2,:,:),[np nElems]);
    coef_mass = (ke.^2).*ccge.*param1.*param2;
    coef_mass = bsxfun(@times,coef_mass,detJ);
    Me = AMS_I_K * coef_mass;
    Me = reshape(Me,[np np nElems]);
    
    %Diffusion matrix
    coefx = ccge.*(param2./param1);
    mult_aux = invdetJ.*invJ11;
    coefx_1 = bsxfun(@times,coefx,mult_aux.*invJ11);
    coefx_2_3 = bsxfun(@times,coefx,mult_aux.*invJ12);
    coefx_4 = bsxfun(@times,coefx,invdetJ.*invJ12.*invJ12);
    
    coefy = ccge.*(param1./param2);
    mult_aux = invdetJ.*invJ21;
    coefy_1 = bsxfun(@times,coefy,mult_aux.*invJ21);
    coefy_2_3 = bsxfun(@times,coefy,mult_aux.*invJ22);
    coefy_4 = bsxfun(@times,coefy,invdetJ.*invJ22.*invJ22);
    
    Ke = AKS_xi_xi*(coefx_1+coefy_1) + AKS_sym*(coefx_2_3+coefy_2_3) + AKS_eta_eta*(coefx_4+coefy_4);
    Ke = reshape(Ke,[np np nElems]);
    
    %Difussion - mass
    Ke = Me - Ke;
    
    %Source term
    fe = zeros(np,1,nElems);

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

pos0 = find(IF == 0);
IF(pos0) = pos0; %Trick to add zero values. This match dimensions in nefem computations
f = accumarray(IF,f);
nodesInterior = unique(T);
ftotal(nodesInterior) = f(nodesInterior);
K = sparse(I,J,Kreal + sqrt(-1)*Kimag,nOfNodes,nOfNodes);











