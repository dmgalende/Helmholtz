clearvars
close all
setpath()


%% INPUT DATA

% Mesh
filemesh = {
    'Farfield_expansions/data/centered_R4/meshes/mesh_P2_10ppw.dcm'
    'Farfield_expansions/data/centered_R4/meshes/mesh_P2_20ppw.dcm'
    'Farfield_expansions/data/centered_R4/meshes/mesh_P2_40ppw.dcm'
    'Farfield_expansions/data/centered_R4/meshes/mesh_P2_80ppw.dcm' 
    };

% Solution file
solfile = {
    'Farfield_expansions/data/centered_R4/sol_KDFE/KDFE_neum_R4_2pi_P2_10ppw.mat'
    'Farfield_expansions/data/centered_R4/sol_KDFE/KDFE_neum_R4_2pi_P2_20ppw.mat'
    'Farfield_expansions/data/centered_R4/sol_KDFE/KDFE_neum_R4_2pi_P2_40ppw.mat'
    'Farfield_expansions/data/centered_R4/sol_KDFE/KDFE_neum_R4_2pi_P2_80ppw.mat'
    };

% Geometry
rad_int = 1;
rad_ext = 4;
x_ext_offset = 0;
x_int_offset = 0;

% Problem
nwaves = -1;
wavenumber = 2*pi;

% Points per wavelength, only applies to compute the element size
% (otherwise use convergence wrt number of dofs)
ppw = [10, 20, 40, 80];


%% CONVERGENCE

% Some postprocessed variables
size_domain = rad_ext - rad_int;
if nwaves > 0
    wavenumber = 2 * pi * nwaves / size_domain;
else
    nwaves = wavenumber * size_domain / (2 * pi);
end

npoints = length(filemesh);
err = zeros(npoints, 1);
hvec = zeros(npoints, 1);

for i = 1:npoints
    
    %% LOAD MESH AND USEFUL VARIABLES
    
    disp(['File ' num2str(filemesh{i})])
    
    % Load mesh
    filemesh_ = GenerateMatFileFromEZ4U(filemesh{i});
    mesh = load(filemesh_{1});
    delete(filemesh_{1})
    deg = size(mesh.Tb_ext, 2) - 1;
    dof_X = size(mesh.X, 1);
    nodes_rad_int = unique(mesh.Tb_int);
    nodes_rad_ext = unique(mesh.Tb_ext);
    dof_ext = length(nodes_rad_ext);
    
    % Vector of element sizes
    %hvec(i) = size_domain * deg / (nwaves * ppw(i));
    hvec(i) = size(mesh.X, 1);

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
    
    %% COMPUTE L2 ERROR
    
    % Solution file
    load(solfile{i})
    
    % Error computation
    err(i) = computeL2NormBoundary(mesh.refelem, mesh.X, mesh.Tb_ext, u_all(:,end),...
                                   @exactSolScatteringCercle, wavenumber, rad_int);
    %err(i) = computeL2NormBoundary(mesh.refelem, mesh.X, mesh.Tb_int, u_all(:,end),...
    %                               @exactSolScatteringCercle, wavenumber, rad_int);

end
err = real(err);

%% PLOT

hold on
figid = plot(hvec, err, '-o'); ax = gca(); ax.XScale = 'log'; ax.YScale = 'log';
for i = 1:npoints-1
    m = log(err(i+1)/err(i)) / log(hvec(i+1)/hvec(i));
    xm = hvec(i);
    ym = err(i);
    text(xm, ym, num2str(round(m, 2)), 'FontSize', 15)
end





