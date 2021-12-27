function [K,f] = berkhoffMatricesAndVolumeVector...
    (X,T,theReferenceElement,u0Vec,kVec,ccgVec,sigmaVec,omega,inPMLelems)

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

%Memory allocation for the interior elements
nOfElementNodes = size(T,2);
aux_ones = ones(1,nOfElementNodes);
allocation = nOfElementNodes^2*nOfElements;
Iint = zeros(allocation,1);                            % allocation??1
Jint = Iint;                                              % allocation??1
Kreal_int = Iint;                                          % allocation??1
Kimag_int = Iint;                                          % allocation??1

%Ordering the elements: PML first
nOfPMLelements = length(inPMLelems);
inINTelems = setdiff(1:nOfElements,inPMLelems);
oElems = [inPMLelems inINTelems];

%Memory allocation for the PML elements
allocation_pml = nOfElementNodes^2*nOfPMLelements;
Ipml = zeros(allocation_pml,1);                              
Jpml = Ipml;                                          
Kreal_pml = Ipml;                                       
Kimag_pml = Ipml;

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
    sigmae = reshape(sigmaVec_t(:,Te),2,nOfElementNodes,1,nElems);  % nD??nOfElementNodes??1??e
    sigmae = permute(sigmae,[2 1 3 4]);                             % nOfElementNodes??nD??1??e
    
    % Assembling parameters
    Te = reshape(Te,1,nOfElementNodes,nElems);
    Te_transp = permute(Te,[2 1 3]);
    aux_row = Te_transp(:,aux_ones,:);
    aux_col = Te(aux_ones,:,:);
    indK = (iniElem-1)*nOfElementNodes^2+1:iElem*nOfElementNodes^2;

    % Elemental Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iElem <= nOfPMLelements %PML elements --> where f~=0 and param > 1
        [k,ccg,param] = computeParameters_PML(N,ccge,ke,sigmae,omega);

        Ke = computeMeMinusKe_tensorProduct_PML(...
            Xe_t,k,param,ccg,Nxieta,np,nd,nD,ngauss,nElems,...
            A_I_K,AM_I_K);
        
        % Assembling for PML area 
        Ipml(indK) = aux_row(:);
        Jpml(indK) = aux_col(:);
        Kreal_pml(indK) = real(Ke);
        Kimag_pml(indK) = imag(Ke);
        
    else %Interior elements --> where f=0 and param=1
        [k,ccg] = computeParameters(N,ccge,ke);

        Ke = computeMeMinusKe_tensorProduct(...
            Xe_t,k,ccg,Nxieta,np,nd,nD,ngauss,nElems,...
            A_I_K,AM_I_K);
        
        % Assembling for interior area 
        Iint(indK) = aux_row(:);
        Jint(indK) = aux_col(:);
        Kreal_int(indK) = real(Ke);
        Kimag_int(indK) = imag(Ke);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update
    iniElem = iElem + 1;
end

pospml = Ipml > 0;
posint = Iint > 0;
Kpml = sparse(Ipml(pospml),Jpml(pospml),Kreal_pml(pospml) + sqrt(-1)*Kimag_pml(pospml),nOfNodes,nOfNodes);
Kint = sparse(Iint(posint),Jint(posint),Kreal_int(posint) + sqrt(-1)*Kimag_int(posint),nOfNodes,nOfNodes);
K = Kpml + Kint;

% Interpolated RHS vector (only inside the PML area).
f = Kpml * u0Vec;


function [k,ccg] = computeParameters(N,ccge,ke)

ccg = sum(bsxfun(@times,N,ccge),1);                         % 1??1??ngauss??e
k = sum(bsxfun(@times,N,ke),1);                             % 1??1??ngauss??e

function [k,ccg,param] = computeParameters_PML(N,ccge,ke,sigmae,omega)

ccg = sum(bsxfun(@times,N,ccge),1);                         % 1??1??ngauss??e
k = sum(bsxfun(@times,N,ke),1);                             % 1??1??ngauss??e
sigma = sum(bsxfun(@times,N,sigmae),1);                     % 1??2??ngauss??e
param = 1 + sqrt(-1)*sigma/omega;                           % 1??2??ngauss??e


function C = computeMeMinusKe_tensorProduct(...
    Xe_t,k,ccg,Nxieta,np,nd,nD,ngauss,nElems,...
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


function C = computeMeMinusKe_tensorProduct_PML(...
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

















