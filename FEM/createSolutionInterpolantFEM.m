function F = createSolutionInterpolantFEM(X,T,u,referenceElement,nDegRef)

% Plotting element (equal spaced points)
nodes = [];
h = 1/nDegRef;
for j = 0:nDegRef
    i = (0:nDegRef-j)';
    aux = j*ones(size(i));
    nodes = [nodes; [i aux]*h];
end
nodes = 2*nodes - 1;
npoints = size(nodes,1);
 
% Vandermonde matrix
nOfElements = size(T,1);
coordRef = referenceElement.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = referenceElement.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');
 
% Compute shape functions at interpolation points
shapeFunctions = zeros(npoints,nOfNodes);
for ipoint = 1:npoints
    p = orthopoly2D(nodes(ipoint,:),nDeg);
    shapeFunctions(ipoint,:) = (invV*p)';
end

% Interpolate solution and position at interpolation points
Xplot = zeros(nOfElements*npoints,2);
uplot = zeros(nOfElements*npoints,1);
fin = 0;
for ielem = 1:nOfElements
    Te = T(ielem,:);
    ind = (fin+1):(fin+npoints);
    Xplot(ind,:) = shapeFunctions*X(Te,:);
    uplot(ind) = shapeFunctions*u(Te);
    fin = fin + npoints;
end

% Set meshfree interpolant (only valid in Matlab 2013)
F = scatteredInterpolant(Xplot(:,1),Xplot(:,2),uplot,'linear','none');



