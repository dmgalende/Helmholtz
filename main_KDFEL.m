clearvars
close all
setpath_FE()


%% INPUT DATA

% Mesh
filemesh = 'Farfield_expansions/data/centered_R4/meshes/mesh_P2_20ppw.dcm';

% Savefile
savefile = 'KDFE_neum_R4_2pi_P2_20ppw.mat';

% Geometry
rad_int = 1;
rad_ext = 4;

% Problem
nwaves = -1;  %define number of waves in the domain, set -1 to use wavenumber value instead
wavenumber = 2*pi;  %only affects if nwaves = -1
inc_angle = 0;

% Boundary condition
bc_type           = 1; %0: soft scattering (Dirichlet BC), 1: partial reflection (Neumann/Robin BC)
bc_reflectionCoef = 0; %only affects if bc_type = 1

% Order of double Karp's single farfield expansion BC
nL = 1:13;

% Which order of family of solutions are going to be depicted
order_to_plot = nL(end);

% Other options
do_renumbering = true;          %use symrcm renumerator before solving the system (recommended)
sort_mesh = false;              %sort the mesh numbering to move the selected boundary nodes at the end
sort_mesh_by_boundary = 'ext';  %only affects if sort_mesh = true


%% STORE DATA

data = struct(...
    'filemesh', filemesh, ...
    'savefile', savefile, ...
    'rad_int', rad_int, ...
    'rad_ext', rad_ext, ...
    'nwaves', nwaves, ...
    'wavenumber', wavenumber, ...
    'inc_angle', inc_angle, ...
    'bc_type', bc_type, ...
    'bc_reflectionCoef', bc_reflectionCoef, ...
    'do_renumbering', do_renumbering, ...
    'sort_mesh', sort_mesh, ...
    'sort_mesh_by_boundary', sort_mesh_by_boundary ...
);


%% LOAD MESH AND USEFUL VARIABLES

% Some postprocessed variables (using unitary wave velocity)
size_domain = rad_ext - rad_int;
if nwaves > 0
    wavenumber = 2 * pi * nwaves / size_domain;
end

% Load mesh
filemesh_ = GenerateMatFileFromEZ4U(filemesh);
mesh = load(filemesh_{1});
delete(filemesh_{1})

% Reference space
deg = size(mesh.Tb_ext, 2) - 1;
mesh.refelem = createReferenceElement(1, (deg+1)*(deg+2)/2, []);

% Rotate exterior boundary elements to have the boundary edge at face 1
mapFace2 = load(['rotTriangle_' num2str(mesh.elemInfo.nOfNodes) 'nodes_face2']);
mapFace3 = load(['rotTriangle_' num2str(mesh.elemInfo.nOfNodes) 'nodes_face3']);
for belem = 1:size(mesh.Tb_ext, 1)
    elem = mesh.elementFaceInfo.ext(belem,1);
    face = mesh.elementFaceInfo.ext(belem,2);
    if face == 2
        mesh.T(elem, :) = mesh.T(elem, mapFace2.rotationMap);
    elseif face == 3
        mesh.T(elem, :) = mesh.T(elem, mapFace3.rotationMap);
    end
    mesh.elementFaceInfo.ext(belem,2) = 1; % new face = 1
end

% Ordering mesh DOFs with selected boundary at the end of the numbering
if sort_mesh
    new_order = orderMesh(mesh, sort_mesh_by_boundary, false);
    mesh.T = new_order(mesh.T);
    mesh.Tb_int = new_order(mesh.Tb_int);
    mesh.Tb_ext = new_order(mesh.Tb_ext);
    mesh.X(new_order, :) = mesh.X;
end

% Variables associated to DOFs
dof_X = size(mesh.X, 1);
nodes_rad_int = unique(mesh.Tb_int);
nodes_rad_ext = unique(mesh.Tb_ext);
dof_ext = length(nodes_rad_ext);

% 2D shape functions on 1D gauss points for face 1
V = Vandermonde_LP(deg, mesh.refelem.NodesCoord);
[L,U,P] = lu(V');
nOfGaussPoints = length(mesh.refelem.IPcoordinates1d);
for i = 1:nOfGaussPoints
    x = [mesh.refelem.IPcoordinates1d(i), -1]; %this is face 1
    p = orthopoly2D_deriv_xieta(x, deg);
    N = U\(L\(P*p));
    mesh.refelem.N_onFace1(i,:) = N(:,1);
end

% 1D angular mesh for farfield functions constructed from boundary mesh (thus it has a periodic connectivity)
mesh.T_ext = restartNumbering(mesh.Tb_ext);
[mesh.X_ext, x_rad] = cart2pol(mesh.X(nodes_rad_ext,1), mesh.X(nodes_rad_ext,2));

% Periodicity will be applied at theta = 0 (thus boundary mesh Tb_ext must have a node!)
periodic_node_pos = find(abs(mesh.X_ext) < 1e-10);
check_val = mesh.X_ext(periodic_node_pos);
if check_val < 0 % be sure that conversion from cartesian to polar provides positive sign!
    mesh.X_ext(periodic_node_pos) = abs(check_val); 
end
periodic_elem = find(mesh.T_ext(:,end) == periodic_node_pos);

% Angular coordinate from [-pi, pi] to [0, 2*pi]
neg_theta = mesh.X_ext < 0;
mesh.X_ext(neg_theta) = 2 * pi + mesh.X_ext(neg_theta);

% Planar incoming wave
wavenumber_vector = [cos(inc_angle); sin(inc_angle)];
uinc = exp( 1i*wavenumber*mesh.X*wavenumber_vector );


%% DISCRETIZATION

% Mass and stiffnes matrices
auxones = ones(dof_X, 1);
[K, M] = PGDberkhoffVolumeMatrices(...
    mesh.X,...
    mesh.T,...
    mesh.refelem,...
    auxones,... %for Kint
    auxones,... %for Kxpml
    auxones,... %for Kypml
    auxones,... %for varMint, varMpml
    []);        %No PML elements
K = K{1};
M = M{1};
    
% Partial/total reflection on physical boundary
if bc_type == 1
    C = berkhoffDampingMatrix(...
        mesh.X,...
        mesh.Tb_int,...
        mesh.refelem,...
        wavenumber * auxones,...
        auxones);
    fc = berkhoffReflectingVector_planarWave(...
         mesh.X,...
         mesh.Tb_int,... 
         mesh.refelem,...
         wavenumber * wavenumber_vector.',...
         1,... %uinc amplitude
         wavenumber * auxones,...
         auxones,...
         bc_reflectionCoef);
else
    C = sparse(dof_X, dof_X);
    fc = zeros(dof_X, 1);
end

% Coupled term on exterior boundary for farfield functions
H = coupledMassMatrix_u(...
    mesh.X,...
    mesh.Tb_ext,...
    mesh.T_ext,...
    mesh.refelem);

% Coupled term on exterior boundary for farfield functions
B = rad_ext * coupledMassMatrix_f(...
    mesh.X_ext,...
    mesh.T_ext,...
    mesh.T,...
    mesh.refelem,...
    periodic_elem,...
    mesh.elementFaceInfo.ext);

% Mass matrix for farfield functions (1D problem)
Mabc = rad_ext * massMatrix1D(mesh.X_ext, mesh.T_ext, mesh.refelem, periodic_elem);

% Stiffness matrix for farfield functions (1D problem)
Kabc = rad_ext * stiffMatrix1D(mesh.X_ext, mesh.T_ext, mesh.refelem, periodic_elem);


%% SYSTEM AND SOLUTION

% Farfield parameters
H0 = besselh(0, wavenumber * rad_ext);
H1 = besselh(1, wavenumber * rad_ext);
[beta, beta2, mu, mu2, delta, delta2] = eval_KDFE_coefs(0, wavenumber, rad_ext, H0, H1);
acoef = delta + (1/rad_ext)*mu + wavenumber^2*beta;
acoef2 = delta2 + (1/rad_ext)*mu2 + wavenumber^2*beta2;
betaR = beta / rad_ext^2;
beta2R = beta2 / rad_ext^2;

% Zero matrices used to dynamically expand the system
zero_vec2 = zeros(2*dof_ext, 1);
zero_mat = sparse(dof_ext, dof_ext);
zero_mat2 = sparse(2*dof_ext, 2*dof_ext);
zero_mat_u = sparse(dof_ext, dof_X);
zero_mat2_u = sparse(2*dof_ext, dof_X);

% Helmholtz system
A = -K + wavenumber^2 * M + 1i*wavenumber*bc_reflectionCoef*C;
if bc_type == 0
    uinc_rad_int = uinc(nodes_rad_int);  
    [A, fc] = applyEssentialBC(A, fc, -uinc_rad_int, nodes_rad_int);
end

% Initial coupled system with the first farfield function
sys_mat = [A , mu * H , mu2 * H
           B , -beta * Mabc, -beta2 * Mabc 
           zero_mat_u, acoef * Mabc - betaR * Kabc, acoef2 * Mabc - beta2R * Kabc];
sys_vec = [fc ; zero_vec2];

% Inits
nnL = length(nL);
u_all = zeros(dof_X, nnL);
f_all = cell(nnL, 1);
ini = 1;
fin = nL(1) - 1;

% Compute solutions for different farfield expansion order
for l = 1:nnL

    disp(['Order = ' num2str(nL(l)) ' with limits ini:fin = ' num2str(ini) ':' num2str(fin)])

    % Update system with the rest of farfield functions (dynamic update, if slow consider allocate)
    for i = ini:fin

        % Coefs for current farfield function
        [beta, beta2, mu, mu2, delta, delta2] = eval_KDFE_coefs(i, wavenumber, rad_ext, H0, H1);
        acoef = delta + (1/rad_ext)*mu + wavenumber^2*beta;
        acoef2 = delta2 + (1/rad_ext)*mu2 + wavenumber^2*beta2;
        betaR = beta / rad_ext^2;
        beta2R = beta2 / rad_ext^2;

        % New column block matrix
        new_col = [mu * H, mu2 * H
                   -beta * Mabc, -beta2 * Mabc
                   acoef * Mabc - betaR * Kabc, acoef2 * Mabc - beta2R * Kabc];
        for j = 1:i-1
            new_col = [new_col ; zero_mat2];
        end

        % New row block matrix
        new_row = zero_mat2_u;
        for j = 1:i-1
            new_row = [new_row , zero_mat2];
        end
        aux_mat = [-(i-1)^2 * Mabc + Kabc, zero_mat ; zero_mat, i^2 * Mabc - Kabc];
        new_row = [new_row , aux_mat];

        % Update system
        aux_mat = [zero_mat, 2*i * Mabc ; 2*i * Mabc, zero_mat];
        sys_mat = [sys_mat , new_col
                   new_row , aux_mat];
        sys_vec = [sys_vec ; zero_vec2];       
    end

    % Update loop limits
    ini = fin + 1;
    if l < nnL, fin = nL(l+1) - 1; end
    
    % Renumbering and FEM solution
    if do_renumbering 
        idx_symrcm = symrcm(sys_mat);
        U = zeros(size(sys_vec));
        U(idx_symrcm) = sys_mat(idx_symrcm, idx_symrcm) \ sys_vec(idx_symrcm);
    else
        U = sys_mat \ sys_vec;
    end
    u_all(:,l) = U(1:dof_X);
    f_all{l} = reshape(U(dof_X+1:end), dof_ext, 2*nL(l));

end


%% SOLUTION PLOTS

u = u_all(:, order_to_plot);
f = f_all{order_to_plot};

disp(['Plotting solution computed with order ', num2str(order_to_plot) , '...'])
figure, plotSolution(mesh.X, mesh.T, real(u), mesh.refelem); title('Scattered wave - Real part');
figure, plotSolution(mesh.X, mesh.T, abs(u), mesh.refelem); title('Scattered wave - Wave height');
figure, plotSolution(mesh.X, mesh.T, real(u+uinc), mesh.refelem); title('Total wave - Real part');
figure, plotSolution(mesh.X, mesh.T, abs(u+uinc), mesh.refelem); title('Total wave - Wave height');

x_theta = cart2pol(mesh.X(nodes_rad_int,1), mesh.X(nodes_rad_int,2));
[x, p] = sort(x_theta);
figure, plot(x, abs(u(nodes_rad_int(p)))), xlabel('Theta'), title('Solution on obstacle (ABS)')

[x, p] = sort(mesh.X_ext);
figure, plot(x, abs(u(nodes_rad_ext(p)))), xlabel('Theta'), title('Solution on artificial boundary (ABS)')
for i = 1:order_to_plot
    f0 = f(:, 2*i-1);
    g0 = f(:, 2*i);
    
    figure
    subplot(1,2,1)
    plot(x, abs(f0(p))), xlabel('Theta'), title(['Farfield function F' num2str(i-1) ' (ABS)'])
    subplot(1,2,2)
    plot(x, abs(g0(p))), xlabel('Theta'), title(['Farfield function G' num2str(i-1) ' (ABS)'])
end


%% SAVE RESULTS AND DATA

if ~isempty(savefile), save(savefile, 'data', 'u_all', 'f_all'); end


