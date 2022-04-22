function [u, f, solverdata] = solve_KDFEL_noDecomposition(mesh, system, uinc, data)

% Input variables
K = system.K;
M = system.M;
H = system.H;
B = system.B;
Mabc = system.Mabc;
Kabc = system.Kabc;
wavenumber = data.problem.wavenumber;
nterms = data.farfield.nterms;
rad_ext = mesh.rad_ext;

% Some inits
dof_X = size(mesh.X, 1);
nodes_rad_int = unique(mesh.Tb_int);
nodes_rad_ext = unique(mesh.Tb_ext);
dof_ext = length(nodes_rad_ext);
fc = zeros(dof_X, 1);

% Farfield parameters
H0 = besselh(0, wavenumber * rad_ext);
H1 = besselh(1, wavenumber * rad_ext);
[beta, beta2, mu, mu2, delta, delta2] = eval_KDFE_coefs(0, wavenumber, rad_ext, H0, H1);
acoef = delta + (1/rad_ext)*mu + wavenumber^2*beta;
acoef2 = delta2 + (1/rad_ext)*mu2 + wavenumber^2*beta2;
betaR = beta / rad_ext^2;
beta2R = beta2 / rad_ext^2;
coef0 = wavenumber * rad_ext;

% Zero matrices used to dynamically expand the system
zero_vec2 = zeros(2*dof_ext, 1);
zero_mat = sparse(dof_ext, dof_ext);
zero_mat2 = sparse(2*dof_ext, 2*dof_ext);
zero_mat_u = sparse(dof_ext, dof_X);
zero_mat2_u = sparse(2*dof_ext, dof_X);

% Helmholtz system
A = -K + wavenumber^2 * M;
uinc_rad_int = uinc(nodes_rad_int);
if nargout == 3
    solverdata = struct();
    [A, fc, solverdata.dirichlet_A_columns] = applyEssentialBC(A, fc, -uinc_rad_int, nodes_rad_int);
else
    [A, fc] = applyEssentialBC(A, fc, -uinc_rad_int, nodes_rad_int);
end

% Initial coupled system with the first farfield function
sys_mat = [A , mu * H , mu2 * H
           B , -beta * Mabc, -beta2 * Mabc 
           zero_mat_u, acoef * Mabc - betaR * Kabc, acoef2 * Mabc - beta2R * Kabc];
sys_vec = [fc ; zero_vec2];

% Update system with the rest of farfield functions (dynamic update, if slow consider allocate)
for i = 1:nterms-1

    % Coefs for current farfield function
    [beta, beta2, mu, mu2, delta, delta2] = eval_KDFE_coefs(i, wavenumber, rad_ext, H0, H1);
    acoef = delta + (1/rad_ext)*mu + wavenumber^2*beta;
    acoef2 = delta2 + (1/rad_ext)*mu2 + wavenumber^2*beta2;
    betaR = beta / rad_ext^2;
    beta2R = beta2 / rad_ext^2;
    coef = coef0 ^ i;

    % New column block matrix
    new_col = coef * [mu * H, mu2 * H
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
    aux_mat = [zero_mat, 2*i*coef0 * Mabc ; 2*i*coef0 * Mabc, zero_mat];
    sys_mat = [sys_mat , new_col
               new_row , aux_mat];
    sys_vec = [sys_vec ; zero_vec2];       
end

% Renumbering and FEM solution
U = zeros(size(sys_vec));
idx = symrcm(sys_mat);
if nargout == 3
    solverdata.idx = idx;
    solverdata.dA = decomposition(sys_mat(idx, idx), 'lu');
    U(idx) = solverdata.dA \ sys_vec(idx);
else
    U(idx) = sys_mat(idx, idx) \ sys_vec(idx);
end   
u = U(1:dof_X);
f = reshape(U(dof_X+1:end), dof_ext, 2*nterms);

% Change of variables to retrieve original functions
coef_l = 0;
for idx = 1:2:2*nterms
    f(:, idx:idx+1) = f(:, idx:idx+1) * coef0^coef_l;
    coef_l = coef_l + 1;
end

    
    
