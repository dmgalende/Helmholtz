home, clear all

useFemNodes = false; %use fem (true) or nefem (false) nodes to compute the boundary
fem         = load('ex_p2_readed.mat');
nefem       = load('ex_p1_nefemMesh.mat');
nDegRef     = 30;

nefemfields = fieldnames(nefem);
boundaries = {};
cont = 1;
for i = 1:length(nefemfields)
    ifield = nefemfields{i};
    if length(ifield) >= 3 && strcmp(ifield(1:3),'Tb_')
        boundaries{cont} = ifield;
        cont = cont + 1;
    end
end

%Plotting element (equal spaced points)
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
referenceElement = createReferenceElement(1,fem.elemInfo.nOfNodes);
coordRef = referenceElement.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = referenceElement.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');

%Compute shape functions at interpolation points
shapeFunctions = zeros(npoints,nOfNodes);
for ipoint = 1:npoints
    p = orthopoly2D(nodes(ipoint,:),nDeg);
    shapeFunctions(ipoint,:) = (invV*p)';
end

if useFemNodes, X = fem.X; else X = nefem.X; end
for i = 1:length(boundaries);
    iname = boundaries{i};
    T = nefem.(iname);
    
    for ielem = 1:size(T,1);
        inodes = 1:nDegRef+1; %nefem element has the curved boundary in the 1st face!
        Te = T(ielem,:);
        Xe = shapeFunctions*X(Te,:);
        Xplot = Xe(inodes,:);
        
        hold on
        plot(Xplot(:,1),Xplot(:,2),'k--','linewidth',2);
    end
end



