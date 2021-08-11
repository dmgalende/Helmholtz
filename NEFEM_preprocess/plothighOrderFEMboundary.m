home, clear all

fem      = load('ex_p15_readed.mat');
nDegRef  = 30;

boundaries = fieldnames(fem.elementFaceInfo);

%Discretize the three faces of the reference element (equal spaced points)
aux = linspace(-1,1,nDegRef+1)';
auxones = ones(nDegRef+1,1);
nodes(:,:,1) = [aux,-auxones];
nodes(:,:,2) = [flipud(aux),-flipud(aux)];
nodes(:,:,3) = [-auxones,flipud(aux)];

% Vandermonde matrix
referenceElement = createReferenceElement(1,fem.elemInfo.nOfNodes);
coordRef = referenceElement.NodesCoord;
nOfNodes = size(coordRef,1);
nDeg = referenceElement.degree;
V = Vandermonde_LP(nDeg,coordRef);
invV = inv(V');

%Compute shape functions at interpolation points
shapeFunctions = zeros(nDegRef+1,nOfNodes,3);
for i = 1:3
    inodes = nodes(:,:,i);
    for ipoint = 1:nDegRef+1
        p = orthopoly2D(inodes(ipoint,:),nDeg);
        shapeFunctions(ipoint,:,i) = (invV*p)';
    end
end

%Interpolate boundary and plot
X = fem.X;
T = fem.T;
for i = 1:length(boundaries)
    iname = ['Tb_' boundaries{i}];
    Tb = fem.(iname);
    ielemFaceInfo = fem.elementFaceInfo.(boundaries{i});
    
    for ielem = 1:size(Tb,1)
        iface = ielemFaceInfo(ielem,2);
        elem2d = ielemFaceInfo(ielem,1);  
        Te = T(elem2d,:);
        Xe = shapeFunctions(:,:,iface)*X(Te,:);
        
        hold on
        plot(Xe(:,1),Xe(:,2),'k--','linewidth',2);
    end
end



