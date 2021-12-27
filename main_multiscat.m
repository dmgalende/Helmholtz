clearvars
close all
setpath_FE()


%% DATA

% Obstacle mesh files
data.mesh.files = {
    'Farfield_expansions/data/centered_R2/meshes/mesh7_P4_10ppw.dcm'
    'Farfield_expansions/data/centered_R2/meshes/mesh7_P4_10ppw.dcm'
    %'Farfield_expansions/data/centered_R2/meshes/mesh7_P4_10ppw.dcm'
    %'Farfield_expansions/data/centered_R2/meshes/mesh7_P4_10ppw.dcm'
    };
nobj = length(data.mesh.files);

% External radius of each domain
rad = 2 * ones(nobj, 1);

% Global coordinate center of each domain
dx = 2;
data.mesh.center = zeros(nobj, 2);
for i = 2:nobj
    data.mesh.center(i,1) = data.mesh.center(i-1,1) + rad(i-1) + dx + rad(i);
end

% Problem setup
data.problem.wavenumber = 2 * pi;
data.problem.max_iter = 50;
data.problem.inc_angle = 0;
data.problem.max_error = 1e-10;

% Farfield expansions setup
data.farfield.nterms = 15;

% Compute analytical solution for circle obstacles
analytical.compute = true;
analytical.M = 100;
analytical.Nfar = 25;
analytical.acoustic_hard = 0;


%% MESHES AND SYSTEM MATRICES

mesh = cell(nobj, 1);
K = cell(nobj, 1);
M = cell(nobj, 1);
H = cell(nobj, 1);
B = cell(nobj, 1);
Mabc = cell(nobj, 1);
Kabc = cell(nobj, 1);
for i = 1:nobj
    mesh{i} = read_mesh(data.mesh.files{i});
    [K{i}, M{i}, H{i}, B{i}, Mabc{i}, Kabc{i}] = compute_system(mesh{i});
end
system = struct('K', K, 'M', M, 'H', H, 'B', B, 'Mabc', Mabc, 'Kabc', Kabc);


%% MULTISCATTERING ITERATIVE SCHEME

% Inits
sol_interior = cell(nobj, 1);
sol_farfield = cell(nobj, 1);
sol_incident = cell(nobj, 1);
wavenumber_vector = [cos(data.problem.inc_angle); sin(data.problem.inc_angle)];
for i = 1:nobj
    sol_interior{i} = zeros(size(mesh{i}.X, 1), 1);
    sol_farfield{i} = zeros(size(mesh{i}.X_ext, 1), 2*data.farfield.nterms);
    sol_incident{i} = exp(1i * data.problem.wavenumber * mesh{i}.X * wavenumber_vector);
end
solution = struct('u', sol_interior, 'f', sol_farfield);
solution_prev = solution;
niter = 0;
convergence = false;
err = zeros(nobj, 1);
indicator = zeros(data.problem.max_iter, 1);

% Iteration loop
while ~convergence && (niter < data.problem.max_iter)

    % Compute solution for each obstacle
    for i = 1:nobj
        
        % Assemble contributions from other obstacles into the incident field for the current obstacle
        uinc = sol_incident{i};
        nodes_int = unique(mesh{i}.Tb_int);
        coord_int = mesh{i}.X(nodes_int, :);
        for j = setdiff(1:nobj, i)
            coord_int_local = coord_int + (data.mesh.center(i,:) - data.mesh.center(j,:));
            uj = evalFarfield(solution_prev(j).f, mesh{j}, coord_int_local, data.problem.wavenumber);
            uinc(nodes_int) = uinc(nodes_int) + uj;
        end
        
        % Solve problem with KDFE boundary conditions
        [solution(i).u, solution(i).f] = solve_KDFEL(mesh{i}, system(i), uinc, data);
        
        % Error on the obstacle wrt the previous solution
        u = solution(i).u(nodes_int);
        u_prev = solution_prev(i).u(nodes_int);
        err(i) = norm(u - u_prev) / norm(u);
        
        % Gauss-Seidel update (disable Jacobi update)
        solution_prev(i) = solution(i);
        
    end
    
    % Convergence criteria
    indicator(niter+1) = max(err);
    convergence = indicator(niter+1) <= data.problem.max_error;
    disp(['Iter ', num2str(niter), ': ', num2str(indicator(niter+1))])
    
    % Iter update
    niter = niter + 1;
    
    % Jacobi update (disable Gauss-Seidel update)
    % solution_prev = solution;

end

% After convergence, ensamble contributions from all scattered fields into each solution
for i = 1:nobj
    for j = setdiff(1:nobj, i)
        coord_local = mesh{i}.X + (data.mesh.center(i,:) - data.mesh.center(j,:));
        uj = evalFarfield(solution(j).f, mesh{j}, coord_local, data.problem.wavenumber);
        solution(i).u = solution(i).u + uj;
    end
end

% Plot convergence
figure
semilogy(0:niter-1, indicator(1:niter), '-o')
title('Convergence')
xlabel('Iteration'), ylabel('Error indicator')


%% ANALYTICAL SOLUTION FOR CYLINDER OBSTACLES

if analytical.compute

    % Data required by the analytical function
    all_rad = zeros(nobj, 1);
    Bx = cell(nobj, 1);
    By = cell(nobj, 1);
    opos = cell(nobj, 1);
    theta = cell(nobj, 1);
    for i = 1:nobj
        nodes_ext = unique(mesh{i}.Tb_ext);
        nodes_int = unique(mesh{i}.Tb_int);
        [~, all_rad(i)] = cart2pol(mesh{i}.X(nodes_int(1), 1), mesh{i}.X(nodes_int(1), 2));
        Bx{i} = mesh{i}.X(nodes_ext, 1);
        By{i} = mesh{i}.X(nodes_ext, 2);
        [theta{i}, opos{i}] = sort(cart2pol(Bx{i}, By{i}), 'ascend');
        Bx{i} = Bx{i}(opos{i}) + data.mesh.center(i, 1);
        By{i} = By{i}(opos{i}) + data.mesh.center(i, 2);
    end

    % Compute analytical solution
    u_exact = NcylExact(data.problem.wavenumber, analytical.acoustic_hard, data.problem.inc_angle, all_rad, ... 
                        data.mesh.center, analytical.M, analytical.Nfar, Bx, By);

    % Comparison plots
    for i = 1:nobj
        figure, hold on
        nodes_ext = unique(mesh{i}.Tb_ext);
        u = solution(i).u(nodes_ext(opos{i}));
        plot(theta{i}, abs(u), 'r-')
        plot(theta{i}, abs(u_exact{i}), 'k--')
        legend('Solution', 'Exact')
        title(['Scattered wave obstacle ', num2str(i), ' - wave height on artificial boundary'])
    end

end


%% SOLUTION PLOTS

fig1 = figure; hold on
fig2 = figure; hold on
fig3 = figure; hold on
fig4 = figure; hold on
for iplot = 1:nobj
    X = mesh{iplot}.X + data.mesh.center(iplot,:);
    figure(fig1.Number)
    plotSolution(X, mesh{iplot}.T, real(solution(iplot).u), mesh{iplot}.refelem);
    figure(fig2.Number)
    plotSolution(X, mesh{iplot}.T, abs(solution(iplot).u), mesh{iplot}.refelem);
    figure(fig3.Number)
    plotSolution(X, mesh{iplot}.T, real(solution(iplot).u + sol_incident{iplot}), mesh{iplot}.refelem);
    figure(fig4.Number)
    plotSolution(X, mesh{iplot}.T, abs(solution(iplot).u + sol_incident{iplot}), mesh{iplot}.refelem);
end
figure(fig1.Number), title('Scattered wave - Real part'), box on
figure(fig2.Number), title('Scattered wave - Wave height'), box on
figure(fig3.Number), title('Total wave - Real part'), box on
figure(fig4.Number), title('Total wave - Wave height'), box on

















