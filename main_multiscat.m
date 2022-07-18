clearvars
close all
setpath_FE()


%% DATA

% Obstacle mesh files
data.mesh.files = {
    'Farfield_expansions/data/centered_R2/meshes/R2_2pi_P2_20ppw.dcm'
    'Farfield_expansions/data/centered_R2/meshes/R2_2pi_P2_20ppw.dcm'
    'Farfield_expansions/data/centered_R2/meshes/R2_2pi_P4_10ppw.dcm'
    %'Farfield_expansions/data/centered_R2/meshes/R2_2pi_P2_40ppw.dcm'
    %'Farfield_expansions/data/centered_R2/meshes/R2_2pi_P2_40ppw.dcm'
    %'Farfield_expansions/data/centered_R2/meshes/R2_2pi_P2_40ppw.dcm'
    };
nobj = length(data.mesh.files);

% Global coordinate center of each domain
dx = 1;
rad = 2 * ones(nobj, 1);
data.mesh.center = zeros(nobj, 2);
for i = 2:nobj
    data.mesh.center(i,1) = data.mesh.center(i-1,1) + rad(i-1) + dx + rad(i);
end
% data.mesh.center = [-2.5 1; 2.5 1; -8 -4; -4 -5; 0 -6; 4 -5; 8 -4];

% General setup
data.setup.doPlots = true;
data.setup.save = false;
data.setup.deg = zeros(nobj, 1);
data.setup.ppw = zeros(nobj, 1);
data.setup.algo = 'gs'; % Jacobi or gs

% Problem setup
data.problem.wavenumber = 2 * pi;
data.problem.max_iter = 30;
data.problem.inc_angle = 0 * pi/180;
data.problem.max_error = 1e-12;

% Farfield expansions setup
data.farfield.nterms = 15;

% Compute analytical solution for circle obstacles
analytical.compute = true;
analytical.M = 150;
analytical.Nfar = 25;
analytical.acoustic_hard = 0;

% Read ppw from file name
for i = 1:nobj
    ppw_idx0 = strfind(data.mesh.files{i}, 'ppw');
    ppw_idx1 = strfind(data.mesh.files{i}, '_');
    data.setup.ppw(i) = str2double(data.mesh.files{i}(ppw_idx1(end)+1:ppw_idx0-1));
end


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
    data.setup.deg(i) = mesh{i}.refelem.degree;
    [K{i}, M{i}, H{i}, B{i}, Mabc{i}, Kabc{i}] = compute_system(mesh{i});
end
system = struct('K', K, 'M', M, 'H', H, 'B', B, 'Mabc', Mabc, 'Kabc', Kabc);


%% MULTISCATTERING ITERATIVE SCHEME

% Inits
sol_interior = cell(nobj, 1);
sol_farfield = cell(nobj, 1);
sol_incident = cell(nobj, 1);
wavenumber_vector = data.problem.wavenumber * [cos(data.problem.inc_angle); sin(data.problem.inc_angle)];
for i = 1:nobj
    sol_interior{i} = zeros(size(mesh{i}.X, 1), 1);
    sol_farfield{i} = zeros(size(mesh{i}.X_ext, 1), 2*data.farfield.nterms);
    sol_incident{i} = exp(1i * mesh{i}.X * wavenumber_vector);
end
solution = struct('u', sol_interior, 'f', sol_farfield);
solution_prev = solution;
niter = 0;
convergence = false;
err = zeros(nobj, 1);
indicator = zeros(data.problem.max_iter, 1);
ddata = cell(nobj, 1);
ddata(:) = {[]};

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
        [solution(i).u, solution(i).f, ddata{i}] = solve_KDFEL(mesh{i}, system(i), uinc, data, ddata{i});
        
        % Error on the obstacle wrt the previous solution
        u = solution(i).u(nodes_int);
        u_prev = solution_prev(i).u(nodes_int);
        err(i) = norm(u - u_prev) / norm(u);
        
        % Gauss-Seidel update
        if strcmpi(data.setup.algo, 'gs')
            solution_prev(i) = solution(i);
        end
        
    end
    
    % Convergence criteria
    indicator(niter+1) = max(err);
    convergence = indicator(niter+1) <= data.problem.max_error;
    disp(['Iter ', num2str(niter), ': ', num2str(indicator(niter+1))])
    
    % Iter update
    niter = niter + 1;
    
    % Jacobi update
    if strcmpi(data.setup.algo, 'jacobi')
        solution_prev = solution;
    end

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
    err = zeros(nobj, 1);
    for i = 1:nobj
        figure, hold on
        nodes_ext = unique(mesh{i}.Tb_ext);
        u = solution(i).u(nodes_ext(opos{i}));
        plot(theta{i}, abs(u), 'r-', 'linewidth', 2.5)
        plot(theta{i}, abs(u_exact{i}), 'k--', 'linewidth', 2.5)
        legend('Solution', 'Exact')
        title(['Scattered wave obstacle ', num2str(i), ' - wave height on artificial boundary'])
        
        % Error computation
        U = zeros(size(mesh{i}.X, 1), 1);
        U0 = zeros(size(mesh{i}.X, 1), 1);
        auxones = ones(size(mesh{i}.X, 1), 1);
        Mext = berkhoffDampingMatrix(...
                mesh{i}.X,...
                mesh{i}.Tb_ext,...
                mesh{i}.refelem,...
                auxones,...
                auxones);
        U(nodes_ext(opos{i})) = u;
        U0(nodes_ext(opos{i})) = u_exact{i};
        e = U - U0;
        err(i) = sqrt(real((e' * Mext * e) / (U0' * Mext * U0)));
        %err(i) = norm(e) / norm(U0); % alternative way with discrete L2 norm
    end
    
    % Save convergece data
    if data.setup.save
        save(['multiscat_errors_n', num2str(nobj), '_k', num2str(data.problem.wavenumber), '_P', ...
             num2str(data.setup.deg(1)), '_', num2str(data.setup.ppw(1)), 'ppw.mat'], 'data', 'err')
    end

end


%% SOLUTION PLOTS

if data.setup.doPlots
    
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
    
end

















