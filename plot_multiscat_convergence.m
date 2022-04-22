clearvars
close all
home

% Load error files
files = {
    %'Farfield_expansions/data/centered_R2/line_objects_dx1/multiscat_errors_n3_k6.2832_P4_10ppw.mat'
    %'Farfield_expansions/data/centered_R2/line_objects_dx1/multiscat_errors_n3_k6.2832_P4_20ppw.mat'
    %'Farfield_expansions/data/centered_R2/line_objects_dx1/multiscat_errors_n3_k6.2832_P4_40ppw.mat'
    %'Farfield_expansions/data/centered_R2/line_objects_dx1/multiscat_errors_n3_k6.2832_P4_80ppw.mat'
    'multiscat_errors_n3_k6.2832_P2_10ppw.mat'
    'multiscat_errors_n3_k6.2832_P2_20ppw.mat'
    'multiscat_errors_n3_k6.2832_P2_40ppw.mat'
    'multiscat_errors_n3_k6.2832_P2_80ppw.mat'
    };

% Object index
object_ind = 1;

% Compute error and element size vector
npoints = length(files);
e = zeros(npoints, 1);
h = zeros(npoints, 1);
for i = 1:npoints
    d = load(files{i});
    e(i) = d.err(object_ind);
    h(i) = 2 * pi * d.data.setup.deg(object_ind) / (d.data.problem.wavenumber * d.data.setup.ppw(object_ind));
end

% Plot
hold on
plotHandle = plot(h, e, 'o-', 'linewidth', 3, 'markersize', 10);
set(gca, 'xScale', 'log')
set(gca, 'yScale', 'log')
set(gca, 'fontsize', 20)
xlabel('Element size')
ylabel('L2 error')
grid on
box on


for i = 1:npoints-1
    m = log(e(i+1)/e(i)) / log(h(i+1)/h(i));
    xm = h(i);
    ym = e(i);
    text(xm, ym, num2str(round(m, 2)), 'FontSize', 15)
end