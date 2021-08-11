function [Ktot,ftot] = applyEssentialBC(K,f,alpha,nodes)

% Input:
%  K: global matrix of the linear system without dirichlet BC
%  f: global vector of the linear system without dirichlet BC
%  alpha: vector of prescribed values
%  nodes: list of the global nodes where the prescribed values will be 
%         assigned
% Output:
%  Ktot: global matrix inclouding prescribed values
%  ftot: global vector inclouding prescribed values
%
% NOTE: the output variables have the same global numbering that K and f

nOfPrescribedNodes = length(nodes);
Ktot = K;
ftot = f;
for inode = 1:nOfPrescribedNodes
    index = nodes(inode);
    value = alpha(inode);
    column = Ktot(:,index);
    column(index) = 0;
    Ktot(index,:) = 0;
    Ktot(:,index) = 0;
    Ktot(index,index) = 1;
    ftot(index) = value;
    ftot = ftot - value*column;
end