function [error eg domArea] = calculateL2ErrorTwoSolutionDifferentP(data_1,data_2,option_solution,ruled)

% calculate the L2 error element by element for two meshes with different p

if data_1.mesh.referenceElement.degree >= data_2.mesh.referenceElement.degree
    data_P = data_1;
    data_p = data_2;
else
    data_P = data_2;
    data_p = data_1;
end
clear data_1 data_2

% incident potential
if isempty(data_p.ip.waveNumberValue)
    if data_p.ip.waveNumberBoundary
        nameConBoundary = data_p.mesh.fieldNames{data_p.mesh.boundaryIndex(...
            data_p.ip.waveNumberBoundary)};
        nodesBoundary = unique(data_p.mesh.(nameConBoundary));
        data_p.ip.waveNumberValue = data_p.bottom.waveNumber(nodesBoundary(1));
    end
end
angle = data_p.ip.direction*(pi/180);
kvector = data_p.ip.waveNumberValue*[cos(angle) sin(angle)];

% Mesh info
T_P = data_P.mesh.T;
X_P = data_P.mesh.X;
referenceElement_P = data_P.mesh.referenceElement;
X_p = data_p.mesh.X;
T_p = data_p.mesh.T;
referenceElement_p = data_p.mesh.referenceElement;
nOfElements = size(T_P,1);
u_p = data_p.solution;
u_P = data_P.solution;

if ruled
    % ruled connectivity
    rul = createRuledConnectivity(X_P,T_P);
else
    % ruled connectivity
    rul = false(nOfElements,1);
    % out of PML elements
    OutOfPMLelements = data_p.OutOfPMLelements;
    rul(OutOfPMLelements) = true;
end
elements = 1:nOfElements;
elements2compute = elements(rul);       % to use the ruled connectivity

% number of nodes per element
nOfNodes_p = size(referenceElement_p.NodesCoord,1); % N of nodes lower mesh
nOfNodes_P = size(referenceElement_P.NodesCoord,1); % N of nodes higher mesh

% Compute shape functions at Gauss points of the higher mesh
coordGauss = referenceElement_P.IPcoordinates; % Gauss points higher mesh
coordRef = referenceElement_p.NodesCoord; % Nodes on the lower mesh
nOfGauss = size(coordGauss,1);
nDeg = referenceElement_p.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');
shapeFunctions = zeros(nOfGauss,nOfNodes_p);
for ipoint = 1:nOfGauss
    [p,dp_dxi,dp_deta] = orthopoly2D_deriv_xieta(coordGauss(ipoint,:),nDeg);
    N = (invV*[p,dp_dxi,dp_deta])';
    shapeFunctions(ipoint,:,1) = N(1,:);
    shapeFunctions(ipoint,:,2) = N(2,:);
    shapeFunctions(ipoint,:,3) = N(3,:);
end

% computation
computation_P = data_P.computation;
computation_p = data_p.computation;

%free memory
clear data_P data_p
loop_ind = 0;

% allocate vectors
error = zeros(nOfElements,1);
solNorm = zeros(nOfElements,1);
domArea = ones(nOfElements,1);
if any(strcmp(computation_p,{'CDG' 'DG'})) % Discontinuous solution
    if any(strcmp(computation_P,{'CDG' 'DG'})) % Discontinuous reference solution
        % Loop in elements
        for ielem = elements2compute
            loop_ind = loop_ind+1;
            Te_P = T_P(ielem,:);
            Xe_P = X_P(Te_P,:);
            Te_p = T_p(ielem,:);
            Xe_p = X_p(Te_p,:);
            ind_p =  (ielem-1)*nOfNodes_p+1:ielem*nOfNodes_p;
            ue_p = u_p(ind_p);
            ind_P =  (ielem-1)*nOfNodes_P+1:ielem*nOfNodes_P;
            ue_P = u_P(ind_P);
            [elemError elemSolNorm elemArea] = calculateElementalError...
                (ue_P,ue_p,referenceElement_P,Xe_P,Xe_p,shapeFunctions,...
                kvector,option_solution,coordGauss,coordRef,nDeg);
            domArea(ielem) = elemArea;
            solNorm(ielem) = elemSolNorm;
            error(ielem) = elemError;
        end
    elseif strcmp(computation_P,'FEM') % Continuous reference solution
        % Loop in elements
        for ielem = elements2compute
            loop_ind = loop_ind+1;
            Te_P = T_P(ielem,:);
            Xe_P = X_P(Te_P,:);
            Te_p = T_p(ielem,:);
            Xe_p = X_p(Te_p,:);
            ind_p =  (ielem-1)*nOfNodes_p+1:ielem*nOfNodes_p;
            ue_p = u_p(ind_p);
            ue_P = u_P(Te_P);
            [elemError elemSolNorm elemArea] = calculateElementalError...
                (ue_P,ue_p,referenceElement_P,Xe_P,Xe_p,shapeFunctions,...
                kvector,option_solution,coordGauss,coordRef,nDeg);
            domArea(ielem) = elemArea;
            solNorm(ielem) = elemSolNorm;
            error(ielem) = elemError;
        end
    end
elseif strcmp(computation_p,{'FEM'}) % Continuous solution
    if any(strcmp(computation_P,{'CDG' 'DG'})) % Discontinuous reference solution
        % Loop in elements
        for ielem = elements2compute
            loop_ind = loop_ind+1;
            Te_P = T_P(ielem,:);
            Xe_P = X_P(Te_P,:);
            Te_p = T_p(ielem,:);
            Xe_p = X_p(Te_p,:);
            ue_p = u_p(Te_p);
            ind_P =  (ielem-1)*nOfNodes_P+1:ielem*nOfNodes_P;
            ue_P = u_P(ind_P);
            [elemError elemSolNorm elemArea] = calculateElementalError...
                (ue_P,ue_p,referenceElement_P,Xe_P,Xe_p,shapeFunctions,...
                kvector,option_solution,coordGauss,coordRef,nDeg);
            domArea(ielem) = elemArea;
            solNorm(ielem) = elemSolNorm;
            error(ielem) = elemError;
        end
    elseif strcmp(computation_P,'FEM') % Continuous reference solution
        % Loop in elements
        for ielem = elements2compute
            loop_ind = loop_ind+1;
            Te_P = T_P(ielem,:);
            Xe_P = X_P(Te_P,:);
            Te_p = T_p(ielem,:);
            Xe_p = X_p(Te_p,:);
            ue_p = u_p(Te_p);
            ue_P = u_P(Te_P);
            [elemError elemSolNorm elemArea] = calculateElementalError...
                (ue_P,ue_p,referenceElement_P,Xe_P,Xe_p,shapeFunctions,...
                kvector,option_solution,coordGauss,coordRef,nDeg);
            domArea(ielem) = elemArea;
            solNorm(ielem) = elemSolNorm;
            error(ielem) = elemError;
        end
    end
end

aux_ind = error~=0;
if strcmpi(option_solution,'FA')
    eg = sqrt(sum(error(aux_ind))/sum(domArea(aux_ind)));
elseif strcmpi(option_solution,'SOL')
    eg = sqrt(sum(error(aux_ind))/sum(solNorm(aux_ind)));
else
    eg = sqrt(sum(error(aux_ind))/sum(domArea(aux_ind)));
end

function [elemError solNorm elemArea] = calculateElementalError...
    (Ue,ue,referenceElement,Xe_P,Xe_p,n,kvector,optSol,coordGauss,coordRef,nDeg)

%Information of the reference element
IPw = referenceElement.IPweights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;

xyg_xi = Nxi*Xe_P;
xyg_eta = Neta*Xe_P;
detJ = xyg_xi(:,1).*xyg_eta(:,2) - xyg_xi(:,2).*xyg_eta(:,1);
dvolu = detJ.*IPw;

% check gauss point position
xyg_P = N*Xe_P;
xyg_p = n(:,:,1)*Xe_p;
xyg_err = (xyg_P - xyg_p);
if max(abs(xyg_err(:)))>1e-10
    gxy0 = transpose(coordGauss);
    [gxy,n,convergence] = recomputeGaussPoint(Xe_p,xyg_P,nDeg,coordRef,gxy0,n);
    if ~convergence
        disp('not converging')
    end
end
Ueg = N*Ue;
ueg = n(:,:,1)*ue;
if strcmpi(optSol,'FA')
    ip_gauss = transpose(exp(sqrt(-1)*kvector*xyg_P'));
    Ueg = abs(Ueg+ip_gauss);
    ueg = abs(ueg+ip_gauss);
end
err_gauss = Ueg-ueg;
elemError = abs(err_gauss'*(dvolu.*err_gauss));
elemArea = sum(dvolu);
solNorm = abs(Ueg'*(dvolu.*Ueg));

function [gxy,n,convergence] = recomputeGaussPoint(Xe_M,xyg_g,nDeg,coordRef,gxy0,n)
ngauss = size(n,1);
res_x = 1;
res_fun = 1;
gxy = gxy0;
tol = 1e-10;
max_iter = 10;
convergence = true;
iter = 1;
while max(abs(res_fun(:)))>tol && max(abs(res_x(:)))>tol && iter<max_iter
    ac = n(:,:,2)*Xe_M; % ng x 2
    bd = n(:,:,3)*Xe_M;
    invDetJ = reshape(1./(ac(:,1).*bd(:,2)-bd(:,1).*ac(:,2)),[1 1 ngauss]);
    ac = permute(ac,[2 3 1]);
    bd = permute(bd,[2 3 1]);
    invJ = bsxfun(@times,([bd(2,1,:) -bd(1,1,:); -ac(2,1,:) ac(1,1,:) ]),invDetJ); % 2 x 2 x ng
    aux = reshape(permute(xyg_g-n(:,:,1)*Xe_M,[2 1]),[ 1 2 ngauss]); % 1 x 2 x ng
    gxy = gxy0 + reshape(sum(bsxfun(@times,invJ,aux),2),[2 ngauss]);
    V = Vandermonde_LP(nDeg,coordRef);
    invV = inv(V');
    [p,dp_dxi,dp_deta] = orthopoly2D_deriv_xieta(gxy',nDeg);
    N = (invV*[p,dp_dxi,dp_deta])';
    n(:,:,1) = N(1:ngauss,:);
    n(:,:,2) = N(ngauss+1:2*ngauss,:);
    n(:,:,3) = N(2*ngauss+1:end,:);
    res_fun = xyg_g - n(:,:,1)*Xe_M;
    res_x = gxy - gxy0;
    gxy0 = gxy;
    iter = iter +1;
end
if iter == max_iter
    convergence = false;
end
