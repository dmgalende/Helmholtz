%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     BERKHOFF EQUATION - HAND FILE                       %
%                     -----------------------------                       %
%                                                                         %
% This script executes a given mesh, bottom, incident potential and       %
% boundary conditions stored in a .mat file created from Berkhoff GUI and %
% computes the solution of Berkhoff equation with them. This solution will%
% be stored into the same path where the data is making use of the same   %
% name but the ending string "_solution" or "_solutionData".              %
%                                                                         %
% This file has been thought to compute the solution of Berkhoff equation %
% without Berkhoff GUI (for instance on an external server). Otherwise    % 
% using the GUI is recommended.                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clf, close all

%-------------------------------- DATA TO BE MODIFIED BY THE USER
% DATA PATH
dataFile = { %Here the data files names ending with .mat    Here the path from current directory ending with /    Here "data" to store complete data or "solution" to store only the solution
             'data_0.5_0.mat'    '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_0.5_1.mat'    '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_0.25_0.mat'   '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_0.25_1.mat'   '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_0.75_0.mat'   '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_0.75_1.mat'   '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_0_0.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_0_1.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_1_0.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_1_1.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_2_0.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_2_1.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_3_0.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_3_1.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_4_0.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             'data_4_1.mat'      '../Harbors_PML/Midscat problem/data_8/'        'data'
             };

%-------------------------------- COMMANDS (DONT MODIFY)
addpath(genpath([pwd '/Callbacks']))
addpath(genpath([pwd '/FEM']))
addpath(genpath([pwd '/NEFEM']))
addpath(genpath([pwd '/CDG']))
addpath(genpath([pwd '/DG']))
for ifile = 1:size(dataFile,1)
    clear data
    path2Open = dataFile{ifile,2};
    file2Open = dataFile{ifile,1};
    what2store = dataFile{ifile,3};
    load([path2Open file2Open])
    computation = data.computation;
    disp(['File: ' file2Open])
    
data.ip.direction = 315; %%%%%%%%%%%%%%

    if strcmp(computation,'FEM')
        data = run_domain(data,[]);
        data = run_boundary(data,[]);
        staticCondensation = data.staticCondensation;
        data2save = {'staticCondensation'};
    elseif strcmp(computation,'NEFEM')
        data = run_domain_NEFEM(data,[]);
        data = run_boundary_NEFEM(data,[]);
        data2save = {};
    elseif strcmp(computation,'CDG')
        data = run_domain_CDG(data,[]);
        data = run_boundary_CDG(data,[]);
        infoFaces = data.infoFaces;
        data2save = {'infoFaces'};
    elseif strcmp(computation,'DG')
        data = preprocess_DG(data,[]);
        data = run_simulation_DG(data);
        infoFaces = data.infoFaces;
        data2save = {'infoFaces'};
    end
    if strcmpi(what2store,'solution')
        sol = data.solution;
        cpuTimeInfo = data.cputime;
        data2save = {data2save{:} 'sol' 'cpuTimeInfo' 'computation'};
        string2save = '_solution.mat';
    else
        data2save = {'data'};
%         string2save = '_solutionData.mat';
string2save = '.mat'; %%%%%%%%%%%%%
    end
data.mesh = rmfield(data.mesh,{'MminusK','fvolume'}); %%%%%%%%%%%%%%%
    save([path2Open file2Open(1:end-4) string2save],data2save{:})
end
    
    
