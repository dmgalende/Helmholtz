function data = run_boundary_NEFEM(data,handle)

% Some useful auxiliar variables and parameters
degree = data.mesh.referenceElement.degree;
angle = data.ip.direction*(pi/180);
kvector = data.ip.waveNumberValue*[cos(angle) sin(angle)];
nOfBoundaries = numel(data.mesh.boundaryNames);
boundaryNodes = zeros(nOfBoundaries,1);
for iboundary = 1:nOfBoundaries
    index = data.mesh.boundaryIndex(iboundary);
    inameCon = data.mesh.fieldNames{index};
    boundaryNodes(iboundary) = length(unique(data.mesh.(inameCon)));
end

% Useful variables to estimate de number of non zero entries in matrices
if data.mesh.elemInfo.type == 1
    valence_nefem = 1;
    vertexPos = 1:3;
elseif data.mesh.elemInfo.type == 0
    valence_nefem = 2;
    vertexPos = 1:4;
end
nOfElementNodes = size(data.mesh.referenceElement.NodesCoord,1);
NNZ_innerNodes = nOfElementNodes;
NNZ_edgeNodes = 2*nOfElementNodes;
NNZ_vertexNodes_nefem = valence_nefem*(nOfElementNodes - degree - 1) + 1;

% Gauss quadrature options
global nIPNurbsEdge
nIPNurbsEdge = 2*degree + 2;

% NRB matrix, damping matrix and reflecting vector from Berkhoff equation (NEFEM)
nOfNodes = size(data.mesh.X,1);
boundaryMatrixSum = spalloc(nOfNodes,nOfNodes,sum(boundaryNodes)*nOfElementNodes); %Poor estimated NNZ entries, not important
boundaryVectorSum = zeros(nOfNodes,1);
t = cputime;
tic
for iboundary = 1:nOfBoundaries
    iname = data.mesh.boundaryNames{iboundary};
    icond = data.BC.values(iboundary);
    iparam = data.BC.parameters{iboundary}{icond};
    inameCon = data.mesh.fieldNames{data.mesh.boundaryIndex(iboundary)};
    itrimInfo = data.mesh.trimmedInfo.(inameCon);
    
    % Number of non zero entries in boundary matrices (NEFEM)
    nOfNodes = boundaryNodes(iboundary);
    nOfElements = size(data.mesh.(inameCon),1);
    nOfInnerNodes = length(data.mesh.referenceElement.innerNodes)*nOfElements;
    nOfLinearNodes = length(unique(data.mesh.(inameCon)(:,vertexPos)));
    nOfEdgeNodes = nOfNodes - nOfLinearNodes - nOfInnerNodes;
    NNZ_nefem = NNZ_vertexNodes_nefem*nOfLinearNodes + NNZ_innerNodes*nOfInnerNodes + NNZ_edgeNodes*nOfEdgeNodes;
    setOutput({['Computation on boundary "' iname '":']},handle)
    
    if icond == 1
        % Do nothing, it is a natural boundary condition of FEM
        setOutput({'       Natural boundary condition'},handle)
    elseif icond == 2
        if iparam
            setOutput({'       Damping matrix...'},handle)
            boundaryMatrixSum = boundaryMatrixSum + sqrt(-1)*iparam*...
                berkhoffDampingMatrix_NEFEM(data.mesh.X,...
                                            data.mesh.(inameCon),...
                                            data.mesh.referenceElement,...
                                            data.mesh.nurbs,...
                                            itrimInfo,...
                                            data.bottom.waveNumber,...
                                            data.bottom.ccg,...
                                            NNZ_nefem);
        end
        setOutput({'       Reflecting vector...'},handle)
        boundaryVectorSum = boundaryVectorSum +...
            berkhoffReflectingVector_NEFEM(data.mesh.X,...
                                           data.mesh.(inameCon),...
                                           data.mesh.referenceElement,...
                                           data.mesh.nurbs,...
                                           itrimInfo,...
                                           kvector,...
                                           data.ip.amplitude,...
                                           data.bottom.waveNumber,...
                                           data.bottom.ccg,...
                                           iparam);
    elseif icond == 3
        setOutput({'       Non reflecting boundary matrix...'},handle)
        boundaryMatrixSum = boundaryMatrixSum +...
            berkhoffNRBCMatrix_NEFEM(data.mesh.X,...
                               data.mesh.(inameCon),...
                               data.mesh.referenceElement,...
                               data.mesh.nurbs,...
                               itrimInfo,...
                               data.bottom.waveNumber,...
                               data.bottom.ccg,...
                               iparam,...
                               NNZ_nefem);
    elseif icond == 4
        setOutput({'       Non reflecting boundary matrix...'},handle)
        boundaryMatrixSum = boundaryMatrixSum +...
            berkhoffNRBCMatrix_NEFEM(data.mesh.X,...
                               data.mesh.(inameCon),...
                               data.mesh.referenceElement,...
                               data.mesh.nurbs,...
                               itrimInfo,...
                               data.bottom.waveNumber,...
                               data.bottom.ccg,...
                               1e99,... %Straight radiation boundary
                               NNZ_nefem); 
    end
end
data.cputime.boundaryToc = toc;
data.cputime.boundary = cputime - t;

% Solution to the linear system
setOutput({'Solving the linear system...'},handle)
A = data.mesh.MminusK + boundaryMatrixSum;
f = boundaryVectorSum - data.mesh.fvolume;
t = cputime;
tic
[L,U,P,Q] = lu(A); %Call UMFPACK
data.solution = Q*(U\(L\(P*f)));
data.cputime.linearSystemToc = toc;
data.cputime.linearSystem = cputime - t;



