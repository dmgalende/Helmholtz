clearvars
close all
setpath()


%% INPUT DATA

% Mesh
filemesh = 'Farfield_expansions/data/centered_R2/meshes/mesh7_P4_20ppw.dcm';

% Solution file
solfile = 'Farfield_expansions/data/centered_R2/sol_KDFE/KDFE_soft_R2_pid4_20ppw.mat';

% Geometry
rad_int      = 1;
rad_ext      = 2;
x_ext_offset = 0;
x_int_offset = 0;

% Problem
nwaves = -1;
wavenumber = pi/4;


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
deg = size(mesh.Tb_ext, 2) - 1;
dof_X = size(mesh.X, 1);
nodes_rad_int = unique(mesh.Tb_int);
nodes_rad_ext = unique(mesh.Tb_ext);
dof_ext = length(nodes_rad_ext);

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

% Reference space
mesh.refelem = createReferenceElement(1, (deg+1)*(deg+2)/2, []);


%% CONVERGENCE PLOT

load(solfile)

[t_ext_0, r_ext_0] = cart2pol(mesh.X(nodes_rad_ext,1) + x_ext_offset, mesh.X(nodes_rad_ext,2));

u0 = zeros(dof_X, 1);
[uscex, A0ex] = ExactSolnFarFieldPatt2D(30, r_ext_0, t_ext_0, wavenumber, rad_int);
u0(nodes_rad_ext) = uscex(1, 1:end-1);

n = size(u_all, 2);
err = zeros(n, 1);
auxones = ones(size(mesh.X, 1), 1);

Mext = berkhoffDampingMatrix(...
        mesh.X,...
        mesh.Tb_ext,...
        mesh.refelem,...
        wavenumber * auxones,...
        auxones);

uk = zeros(dof_X, 1);
for i = 1:n
    uk(nodes_rad_ext) = u_all(nodes_rad_ext, i);
    vec = uk - u0;
    err(i) = (vec' * Mext * vec) / (u0' * Mext * u0);
end
err = sqrt(err);

hold on, plot(0:n-1, err, '-o')
ax = gca(); ax.YScale = 'log';
ylabel('L_2 error on artificial boundary')
xlabel('Farfield approximation order')

[t_ordered, p_ordered] = sort(t_ext_0);
for i = 1:n
    figure
    plot(t_ordered, abs(u0(nodes_rad_ext(p_ordered))), 'b-')
    hold on
    plot(t_ordered, abs(u_all(nodes_rad_ext(p_ordered), i)), 'k--')
    title(['Solution on exterior boundary with k = ', num2str(wavenumber), ' and n = ', num2str(i), ' terms'])
    legend('Exact', 'Numerical')
end

