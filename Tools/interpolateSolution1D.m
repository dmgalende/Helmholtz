function interpu = interpolateSolution1D(u,X,T,referenceElement,interpX)

% Projection of a 1D scalar field u computed in a mesh {X,T,referenceElement} 
% onto a mesh {interpX}

%Mesh info
p = referenceElement.degree;
nElems = size(T,1);

%Vandermonde
V = Vandermonde_LP(p,referenceElement.NodesCoord1d);
[L,U,P] = lu(V');

%Loop in elements of the base mesh
interpu = zeros(size(interpX, 1), size(u, 2));
prepos = 1;
pos = 1;
maxpos = length(interpX);
i = 1;
while i <= nElems && pos <= maxpos
    
    %Search for those interpolated nodes which are included in the element
    x2 = X(T(i,p+1));
    while pos <= maxpos && interpX(pos) <= x2, pos = pos + 1; end
    interpnodes = interpX(prepos:pos-1);
    nOfinterpNodes = length(interpnodes);
    
    %Interpolation nodes in the referenceElement
    xv = X(T(i,[1 p+1]));
    interpnodes_xi = (2*interpnodes - (xv(1)+xv(2))) / (xv(2)-xv(1));
    
    %Compute shape functions at interpolated nodes
    shapeFunctions = zeros(nOfinterpNodes,p+1);
    for j = 1:nOfinterpNodes
        pol = orthopoly1D(interpnodes_xi(j),p);
        shapeFunctions(j,:) = (U\(L\(P*pol)))';
    end
    
    %Interpolate solution at interpolated nodes
    interpu(prepos:pos-1, :) = shapeFunctions*u(T(i,:), :);
    
    %Update
    prepos = pos;
    i = i + 1;
end