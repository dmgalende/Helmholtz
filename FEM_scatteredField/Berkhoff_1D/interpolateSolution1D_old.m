function interpu = interpolateSolution1D(u,X,T,referenceElement,interpX)

% Projection of a 1D scalar field u computed in a mesh {X,T,referenceElement} 
% onto a mesh {interpX}

%Mesh info
p = referenceElement.degree;
nElems = size(T,1);

%Loop in elements of the base mesh
interpu = zeros(size(interpX));
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
    
    %Vandermonde
    V = Vandermonde_LP(p,X(T(i,:)));
    [L,U,P] = lu(V');
    
    %Compute shape functions at interpolated nodes
    shapeFunctions = zeros(nOfinterpNodes,p+1);
    for j = 1:nOfinterpNodes
        pol = orthopoly1D(interpnodes(j),p);
        shapeFunctions(j,:) = (U\(L\(P*pol)))';
    end
    
    %Interpolate solution at interpolated nodes
    interpu(prepos:pos-1) = shapeFunctions*u(T(i,:));
    
    %Update
    prepos = pos;
    i = i + 1;
end