function [nodes, elemsN, elemsS, boundaries, nsubdomains] = importGmsh2D(mesh_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   Import a 2D GMSM mesh file
%
%   Parameters: 
%   In: 
%       mesh_file = file mesh name to be loaded
%   Out:
%       nodes = spatial coordinates of the nodes with dimensions = (number_nodes, 2)
%       elemsN = elements nodal-connectivity with dimensions = (number_elements, 3)
%       elemsS = elements subdomain-connectivity with dimensions = (1,number_elements);
%       nsubdomains = number of subdomains/layers/materials in the mesh
%       boundaries = list of boundary edges with dimensions = (number_boundaries, 3)
%
%   Comments:
%       Each element (row) in boundaries array is composed by three
%       numbers:
%       (id_boundary, node1, node2) 
%
% Author: Octavio Castillo Reyes (octavio.castillo@bsc.es)
% Latest update: July 26th, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----- Import mesh -----%     
% Standard data formats
fstring  = '%s';
finteger = '%d';
% Specific data formats
fnodes   = '%d %f %f %f';
% Open file
Idfile = fopen(mesh_file,'r');
% Read 4 lines (Header)
n = 4; textscan(Idfile,fstring,n,'Delimiter','\n');
% Read number of subdomains
n = 1; input = textscan(Idfile,finteger,n,'Delimiter','\n');
nsubdomains = input{1}(1);
% Read 2 + nsubdomains lines
n = 2 + nsubdomains; textscan(Idfile,fstring,n,'Delimiter','\n');
% Read number of nodes
n = 1; input = textscan(Idfile,finteger,n,'Delimiter','\n');
nNodes = input{1}(1);
% Read spatial positions of the nodes
n = nNodes; input = textscan(Idfile,fnodes,n,'Delimiter','\n');
nodes = (horzcat(input{2},input{3}));
% Read 2 lines
n = 2; input = textscan(Idfile,fstring,n,'Delimiter', '\n');
% Read number of elements
n = 1; input = textscan(Idfile,finteger,n,'Delimiter','\n');
nElems = input{1}(1);
% Read nodal/elements connectivity, edges on boundaries and
% material/elements connectivity
idx_edge = 1;
idx_triangle = 1;
for iline = 1 : nElems
    tline = fgetl(Idfile);
    tmp = str2num(tline);
    if (length(tmp) == 7)   % Edge
        boundaries(idx_edge, 1) = tmp(4);
        boundaries(idx_edge, 2) = tmp(6);
        boundaries(idx_edge, 3) = tmp(7);
        % Update index for edges
        idx_edge = idx_edge + 1;
    elseif (length(tmp) == 8)   % Triangle
        % nodal/elements connectivity
        elemsN(idx_triangle, 1) = tmp(6);
        elemsN(idx_triangle, 2) = tmp(7);
        elemsN(idx_triangle, 3) = tmp(8);
        % material/elements connectivity
        elemsS(idx_triangle) = tmp(4);
        idx_triangle = idx_triangle + 1;
    else
        error('Type of element not supported.');
    end
end
fclose(Idfile);

end