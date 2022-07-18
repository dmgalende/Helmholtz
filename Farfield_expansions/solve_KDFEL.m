function [u, f, solverdata] = solve_KDFEL(mesh, system, uinc, data, ddata)

if (nargin == 5) && ~isempty(ddata)
    [u, f] = solve_KDFEL_withDecomposition(mesh, uinc, data, ddata);
    solverdata = ddata;
else
    [u, f, solverdata] = solve_KDFEL_noDecomposition(mesh, system, uinc, data);
end

