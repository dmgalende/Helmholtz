clearvars
close all
setpath()


%% INPUT DATA

% Mesh
filemesh = 'Farfield_expansions/centered_R2/mesh7_P4_10ppw.dcm';

% Savefile
savefile = 'kk.mat';

% Geometry
rad_int = 1;
rad_ext = 2;

% Problem
nwaves = -1;
wavenumber = 2 * pi;
inc_angle = 0;

% Boundary condition
bc_type           = 1; %0: soft scattering, 1: partial reflection
bc_reflectionCoef = 0; %only affects if bc_type = 1


%% LOAD MESH AND USEFUL VARIABLES

% Some postprocessed variables
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
end

% Reference space
mesh.refelem = createReferenceElement(1, (deg+1)*(deg+2)/2, []);

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

% First order Sommerfeld-like artificial ABC
B = berkhoffNRBCMatrix(...
    mesh.X,...
    mesh.Tb_ext,...
    mesh.refelem,...
    wavenumber * auxones,...
    auxones,...
    rad_ext);

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


%% SYSTEM AND SOLUTION

% System
sys_mat = -K + wavenumber^2 * M + B + 1i*bc_reflectionCoef*C;
sys_vec = fc;
                
% Essential boundary condition
if bc_type == 0
    uinc_rad_int = uinc(nodes_rad_int);  
    [sys_mat, sys_vec] = applyEssentialBC(sys_mat, sys_vec, -uinc_rad_int, nodes_rad_int);
end

% FEM solution
u = sys_mat \ sys_vec;


%% PLOTS

figure, plotSolution(mesh.X, mesh.T, real(u), mesh.refelem); title('Scattered wave - Real part');
figure, plotSolution(mesh.X, mesh.T, abs(u), mesh.refelem); title('Scattered wave - Wave height');
figure, plotSolution(mesh.X, mesh.T, real(u+uinc), mesh.refelem); title('Total wave - Real part');
figure, plotSolution(mesh.X, mesh.T, abs(u+uinc), mesh.refelem); title('Total wave - Wave height');

x_theta = cart2pol(mesh.X(nodes_rad_int,1), mesh.X(nodes_rad_int,2));
[x, p] = sort(x_theta);
figure, plot(x, abs(u(nodes_rad_int(p)))), xlabel('Theta'), title('Solution on obstacle (ABS)')

x_theta = cart2pol(mesh.X(nodes_rad_ext,1), mesh.X(nodes_rad_ext,2));
[x, p] = sort(x_theta);
figure, plot(x, abs(u(nodes_rad_ext(p)))), xlabel('Theta'), title('Solution on artificial boundary (ABS)')


%% Save

if ~isempty(savefile), save(savefile, 'u'); end






