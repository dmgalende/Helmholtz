function u = evalFarfield(f, mesh, x, k)

% Computes the farfield solution u(x) at cartesian points defined in x, by performing the FEM interpolation 
% of the KDFE farfield functions.
% 
% Input args
%   f: matrix of farfield functions of the KDFE boundary conditions.
%   mesh: structure with fields "X_ext" (polar angle of the exterior boundary) and 
%         "refelem" (FEM reference element structure).
%   x: cartesian points where solution will be interpolated.
%   k: wavenumber.
%
% Output arg
%   u: KDFE farfield solution evaluated at cartesian points in x.


% Transform points to polar coordinates in [0, 2*pi] since f depends explicitly on the angle
[theta, rho] = cart2pol(x(:,1), x(:,2));
neg_theta = theta < 0;
theta(neg_theta) = 2 * pi + theta(neg_theta);

% Ascend ordered polar mesh coordinates and f values (required by the interpolation function)
[otheta, opos] = sort(theta, 'ascend');
orho = rho(opos);
[X_ext, opos_ext] = sort(mesh.X_ext, 'ascend');
f_ext = f(opos_ext, :);

% Create unwrapped 1D high-order mesh
X_ext = [X_ext ; 2 * pi];
f_ext = [f_ext ; f_ext(1, :)];
m = (length(X_ext) - 1) / mesh.refelem.degree;
T_ext = create1Dconec(m, mesh.refelem.degree);

% Interpolate farfield functions f at the given angular coordinates using current FEM interpolation
f_otheta = interpolateSolution1D(f_ext, X_ext, T_ext, mesh.refelem, otheta);

% Evaluate farfield
kr = k * orho;
H0 = besselh(0, kr);
H1 = besselh(1, kr);
ou = zeros(size(x,1), 1);
coef = 0;
for i = 1:2:size(f,2)
    ou = ou + (1 ./ kr.^coef) .* (H0 .* f_otheta(:,i) + H1 .* f_otheta(:,i+1));
    coef = coef + 1;
end

% Backward map to original order
u = zeros(size(x,1), 1);
u(opos) = ou;



