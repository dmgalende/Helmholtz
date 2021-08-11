function data = run_domain_NEFEM(data,handle)

% Some useful auxiliar variables and parameters
omega = 2*pi/data.ip.period;
nOfBoundaries = numel(data.mesh.boundaryNames);
degree = data.mesh.referenceElement.degree;

% Celerities and wave number
if ~isfield(data.bottom,'ccg')
    setOutput({'Precalculated celerities not detected'
        'Computing celerities...'},handle)
    if data.constantBottomFlag
        setOutput({'Constant bottom flag is enabled'},handle)
        data.bottom.waveNumber = ones(size(data.mesh.X,1),1)*...
            computeWaveNumber(omega,data.bottom.value(1));
    else
        data.bottom.waveNumber = computeWaveNumber(omega,data.bottom.value);
    end
    data.bottom.ccg = celerities(omega,data.bottom.waveNumber,data.bottom.value);
end

% Definition of initial potential data (constant bottom)
setOutput({'Computing the incident wave...'},handle)
if data.ip.waveNumberBoundary
    nameConBoundary = data.mesh.fieldNames{data.mesh.boundaryIndex(...
        data.ip.waveNumberBoundary)};
    nodesBoundary = getNodes(data.mesh.(nameConBoundary));
    data.ip.waveNumberValue = data.bottom.waveNumber(nodesBoundary(1));
end
angle = data.ip.direction*(pi/180);
kvector = data.ip.waveNumberValue*[cos(angle) sin(angle)];
data.ip.value = data.ip.amplitude * exp(sqrt(-1)*(kvector(1)*data.mesh.X(:,1) ...
    + kvector(2)*data.mesh.X(:,2)));

%Required information for the interior mesh routines
% ne = size(data.mesh.T,1);
% out_e = createOutOfPMLConnectivity(data);
% data.mesh.intT = data.mesh.T(out_e,:);
% data.mesh.inPMLelems = setdiff(1:ne,out_e);
% data.mesh.extT = data.mesh.T(data.mesh.inPMLelems,:);

%Global absorption PML parameter & boundary nodes
data.PMLabsorptionValue = zeros(size(data.mesh.X));
boundaryNodes = zeros(nOfBoundaries,1);
for iboundary = 1:nOfBoundaries
    if strcmp(data.PML{5,iboundary},'on')
        data.PMLabsorptionValue = data.PMLabsorptionValue + ...
            data.PML{1,iboundary};
%         data.PMLabsorptionValue = 0*data.PMLabsorptionValue; %%%%%%%%%%%%%%%%%%%%%%%%
    end
    index = data.mesh.boundaryIndex(iboundary);
    inameCon = data.mesh.fieldNames{index};
    boundaryNodes(iboundary) = length(unique(data.mesh.(inameCon)));
end

% Gauss quadrature options
global nIPNurbsEdge nIPinterior
nIPNurbsEdge = 2*degree + 2;
nIPinterior = 2*degree + 2;

% Useful variables to estimate the number of non zero entries in matrices
if data.mesh.elemInfo.type == 1
    valence_nefem = 1;
    vertexPos = 1:3;
elseif data.mesh.elemInfo.type == 0
    valence_nefem = 2;
    vertexPos = 1:4;
end
nOfElementNodes = size(data.mesh.referenceElement.NodesCoord,1);
NNZ_vertexNodes_nefem = valence_nefem*(nOfElementNodes - degree - 1) + 1;
NNZ_innerNodes = nOfElementNodes;
NNZ_edgeNodes = 2*nOfElementNodes;

% Volume matrices and volume vector from Berkhoff equation (NEFEM)
setOutput({'Computing stiffness, mass matrices and volume vector (NEFEM)...'},handle)
t = cputime;
tic
nOfNodes = size(data.mesh.X,1);
MminusKnefem = spalloc(nOfNodes,nOfNodes,sum(boundaryNodes)*nOfElementNodes); %Poor estimated non zero entries, not important
fnefem = zeros(nOfNodes,1);
nOFBoundaryElements = 0;
for iboundary = 1:nOfBoundaries
    index = data.mesh.boundaryIndex(iboundary);
    inameCon = data.mesh.fieldNames{index};
    itrimInfo = data.mesh.trimmedInfo.(inameCon);
    
    % Number of non zero entries in volume matrices (NEFEM)
    nOfNodes = boundaryNodes(iboundary);
    nOfElements = size(data.mesh.(inameCon),1);
    nOFBoundaryElements = nOFBoundaryElements + nOfElements;
    nOfInnerNodes = length(data.mesh.referenceElement.innerNodes)*nOfElements;
    nOfLinearNodes = length(unique(data.mesh.(inameCon)(:,vertexPos)));
    nOfEdgeNodes = nOfNodes - nOfLinearNodes - nOfInnerNodes;
    NNZ_nefem = NNZ_vertexNodes_nefem*nOfLinearNodes + NNZ_innerNodes*nOfInnerNodes + NNZ_edgeNodes*nOfEdgeNodes;
    
    [MminusKaux,faux] = berkhoffMatricesAndVolumeVector_NEFEM...
        (data.mesh.X,...
         data.mesh.(inameCon),... 
         data.mesh.referenceElement,...
         data.mesh.nurbs,...
         itrimInfo,...
         data.ip.amplitude,...
         kvector,...
         data.bottom.waveNumber,...
         data.bottom.ccg,...
         data.PMLabsorptionValue,...
         omega,...
         NNZ_nefem);
    MminusKnefem = MminusKnefem + MminusKaux;
    fnefem = fnefem + faux;
end

% Volume matrices and volume vector from Berkhoff equation (FEM)
setOutput({'Computing stiffness, mass matrices and volume vector (FEM)...'},handle)
[MminusKfem,ffem] = berkhoffMatricesAndVolumeVector_OLD...
        (data.mesh.X,... 
        data.mesh.T,... %T is connectivity for interior elements
        data.mesh.referenceElement,...
        data.ip.amplitude,...
        kvector,...
        data.bottom.waveNumber,...
        data.bottom.ccg,...
        data.PMLabsorptionValue,...
        omega,...
        nOFBoundaryElements);

data.cputime.volMatToc = toc;
data.cputime.volMat = cputime - t;
data.mesh.MminusK = MminusKfem + MminusKnefem;
data.mesh.fvolume = ffem + fnefem;



