function [u, f] = solve_KDFEL_withDecomposition(mesh, uinc, data, solverdata)

% Input variables
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
coef0 = wavenumber * rad_ext;

% Zero matrices used to dynamically expand the system
zero_vec2 = zeros(2*dof_ext, 1);

% Helmholtz system
uinc_rad_int = uinc(nodes_rad_int);
[~, fc] = applyEssentialBC([], fc, -uinc_rad_int, nodes_rad_int, solverdata.dirichlet_A_columns);

% Initial coupled system with the first farfield function
sys_vec = [fc ; zero_vec2];

% Update system with the rest of farfield functions (dynamic update, if slow consider allocate)
for i = 1:nterms-1
    sys_vec = [sys_vec ; zero_vec2];       
end

% Renumbering and FEM solution
U = zeros(size(sys_vec));
U(solverdata.idx) = solverdata.dA \ sys_vec(solverdata.idx);
u = U(1:dof_X);
f = reshape(U(dof_X+1:end), dof_ext, 2*nterms);

% Change of variables to retrieve original functions
coef_l = 0;
for idx = 1:2:2*nterms
    f(:, idx:idx+1) = f(:, idx:idx+1) * coef0^coef_l;
    coef_l = coef_l + 1;
end

    
    
