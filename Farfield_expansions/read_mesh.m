function mesh = read_mesh(filemesh)

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
nodes_rad_ext = unique(mesh.Tb_ext);
mesh.T_ext = restartNumbering(mesh.Tb_ext);
[mesh.X_ext, rad_ext] = cart2pol(mesh.X(nodes_rad_ext,1), mesh.X(nodes_rad_ext,2));
mesh.rad_ext = rad_ext(1);

% Periodicity will be applied at theta = 0 (thus boundary mesh Tb_ext must have a node!)
periodic_node_pos = find(abs(mesh.X_ext) < 1e-10);
check_val = mesh.X_ext(periodic_node_pos);
if check_val < 0 % be sure that conversion from cartesian to polar provides positive sign!
    mesh.X_ext(periodic_node_pos) = abs(check_val); 
end
mesh.periodic_elem = find(mesh.T_ext(:,end) == periodic_node_pos);

% Angular coordinate from [-pi, pi] to [0, 2*pi]
neg_theta = mesh.X_ext < 0;
mesh.X_ext(neg_theta) = 2 * pi + mesh.X_ext(neg_theta);

