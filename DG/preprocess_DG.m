function data = preprocess_DG(data,handle)

% Some useful auxiliar variables and parameters
nDeg = data.mesh.referenceElement.degree;
omega = 2*pi/data.ip.period;
nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};

% Celerities and wave number
constant_bottom = 0; % 0-variable / 1-constant
if ~isfield(data.bottom,'ccg')
    setOutput({'Precalculated celerities not detected'
        'Computing celerities...'},handle)
    if constant_bottom
        setOutput({'Constant bottom'},handle)
        data.bottom.waveNumber = computeWaveNumber(omega,data.bottom.value(1))*...
            ones(size(data.bottom.value));
        [ccg,c,cg] = celerities(omega,data.bottom.waveNumber(1),data.bottom.value(1));
        data.bottom.ccg = ccg*ones(size(data.bottom.value));
        data.bottom.c = c*ones(size(data.bottom.value));
        data.bottom.cg = cg*ones(size(data.bottom.value));
    else
        data.bottom.waveNumber = computeWaveNumber(omega,data.bottom.value);
        [data.bottom.ccg,data.bottom.c,data.bottom.cg] = celerities(omega,data.bottom.waveNumber,data.bottom.value);
    end
end

% Definition of initial potential data
if data.ip.waveNumberBoundary
    nameConBoundary = data.mesh.fieldNames{data.mesh.boundaryIndex(...
        data.ip.waveNumberBoundary)};
    nodesBoundary = getNodes(data.mesh.(nameConBoundary));
    data.ip.waveNumberValue = data.bottom.waveNumber(nodesBoundary(1));
end

%Global absorption PML parameter
data.PMLabsorptionValue = zeros(size(data.mesh.(nameNodal)));
for iboundary = 1:numel(data.mesh.boundaryNames)
    if strcmp(data.PML{5,iboundary},'on')
        data.PMLabsorptionValue = data.PMLabsorptionValue + ...
            data.PML{1,iboundary};
    end
end

% DG preprocess: infoFaces structure & sector division
if ~isfield(data.mesh,'interiorElements_FreeSpace')
    setOutput({'Mesh preprocess not detected'
        'Performing DG mesh preprocess...'},handle)
    data = mesh_preprocess(data);
end

% Calculating h
Tlinear = data.mesh.(nameCon)(:,1:3);
nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
data.mesh.h = computeMinElementSize(data.mesh.(nameNodal), Tlinear);
% calculating dt
dt = data.mesh.h/(nDeg^2)/max(data.bottom.c); % Runge Kutta
data.sp.dt_stab = dt;

% Inizialization
data = Initialization(data);

% Matrix preprocess: create the matrices for the RE and curved elements
data = matrix_preprocess(data,handle);




