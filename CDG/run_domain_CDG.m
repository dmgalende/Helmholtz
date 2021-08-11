function data = run_domain_CDG(data,handle)

% Some useful auxiliar variables and parameters
omega = 2*pi/data.ip.period;
nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};

% Celerities and wave number
constant_bottom = 1; % 0-variable / 1-constant
if ~isfield(data.bottom,'ccg')
    setOutput({'Precalculated celerities not detected'
        'Computing celerities...'},handle)
    if constant_bottom
        setOutput({'Constant bottom'},handle)
        data.bottom.waveNumber = computeWaveNumber(omega,data.bottom.value(1))*...
            ones(size(data.bottom.value));
        data.bottom.ccg = celerities(omega,data.bottom.waveNumber(1),data.bottom.value(1))*...
            ones(size(data.bottom.value));
    else
        data.bottom.waveNumber = computeWaveNumber(omega,data.bottom.value);
        data.bottom.ccg = celerities(omega,data.bottom.waveNumber,data.bottom.value);
    end
end

% Definition of initial potential data
if data.ip.waveNumberBoundary
    nameConBoundary = data.mesh.fieldNames{data.mesh.boundaryIndex(...
        data.ip.waveNumberBoundary)};
    nodesBoundary = getNodes(data.mesh.(nameConBoundary));
    data.ip.waveNumberValue = data.bottom.waveNumber(nodesBoundary(1));
end
angle = data.ip.direction*(pi/180);
kvector = data.ip.waveNumberValue*[cos(angle) sin(angle)];

%Global absorption PML parameter
data.PMLabsorptionValue = zeros(size(data.mesh.(nameNodal)));
for iboundary = 1:numel(data.mesh.boundaryNames)
    if strcmp(data.PML{5,iboundary},'on')
        data.PMLabsorptionValue = data.PMLabsorptionValue + ...
            data.PML{1,iboundary};
    end
end

% CDG preprocess: infoFaces structure
if ~isfield(data,'infoFaces')
    setOutput({'Precalculated infoFaces not detected'
        'Building CDG infoFaces structure...'},handle)
    data = preprocess_CDG(data);
end

% Calculating h
Tlinear = data.mesh.(nameCon)(:,1:3);
nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
data.mesh.h = computeMinElementSize(data.mesh.(nameNodal), Tlinear);

% P adaptivity preprocess
if isfield(data,'p_adaptivity')
    setOutput({'P adaptivity preprocess (CDG)...'},handle)
    data = p_adaptivity_preprocess(data);
    
    % Volume matrices and volume vector from Berkhoff equation
    setOutput({'Volume computations P-refinement (CDG)...'},handle)
    t = cputime;

    [KminusM,f,R] = berkhoffMatricesAndVolumeVector_CDG_P_adapt...
        (data.mesh.(nameNodal),...
        data.mesh.(nameCon),...
        data.mesh.referenceElement,...
        data.ip.amplitude,...
        kvector,...
        data.bottom.waveNumber,...
        data.bottom.ccg,...
        data.PMLabsorptionValue,...
        omega,data.p_adaptivity);

    data.cputime.volMat_volume = cputime - t;
    t = cputime;

    setOutput({'Interior faces computations P-refinement (CDG)...'},handle)

    H = berkhoffInteriorFaceLoop_CDG_P_adapt...
        (data.mesh.(nameNodal),...
        data.mesh.(nameCon),...
        data.mesh.referenceElement,...
        data.infoFaces,...
        data.bottom.ccg,...
        data.PMLabsorptionValue,...
        omega,R,data.p_adaptivity);
    
else
    % Volume matrices and volume vector from Berkhoff equation
    setOutput({'Volume computations (CDG)...'},handle)
    t = cputime;

    [KminusM,f,R] = berkhoffMatricesAndVolumeVector_CDG...
        (data.mesh.(nameNodal),...
        data.mesh.(nameCon),...
        data.mesh.referenceElement,...
        data.ip.amplitude,...
        kvector,...
        data.bottom.waveNumber,...
        data.bottom.ccg,...
        data.PMLabsorptionValue,...
        omega);

    data.cputime.volMat_volume = cputime - t;
    t = cputime;

    setOutput({'Interior faces computations (CDG)...'},handle)

    H = berkhoffInteriorFaceLoop_CDG...
        (data.mesh.(nameNodal),...
        data.mesh.(nameCon),...
        data.mesh.referenceElement,...
        data.infoFaces,...
        data.bottom.ccg,...
        data.PMLabsorptionValue,...
        omega,R,data.mesh.h);
end

data.cputime.volMat_faces = cputime - t;
data.mesh.KminusMplusH = KminusM + H;
data.mesh.fvolume = f;



