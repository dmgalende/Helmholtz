function data = run_domain(data,handle)

% Some useful auxiliar variables and parameters
omega = 2*pi/data.ip.period;
nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};

% Celerities and wave number
if ~isfield(data.bottom,'ccg')
    setOutput({'Precalculated celerities not detected'
        'Computing celerities...'},handle)
    if data.constantBottomFlag
        setOutput({'Constant bottom flag is enabled'},handle)
        data.bottom.waveNumber = ones(size(data.mesh.(nameNodal),1),1)*...
            computeWaveNumber(omega,data.bottom.value(1));
    else
        data.bottom.waveNumber = computeWaveNumber(omega,data.bottom.value);
    end
    data.bottom.ccg = celerities(omega,data.bottom.waveNumber,data.bottom.value);
end

% Definition of initial potential data
setOutput({'Computing the incident wave...'},handle)
if data.ip.waveNumberBoundary
    nameConBoundary = data.mesh.fieldNames{data.mesh.boundaryIndex(...
        data.ip.waveNumberBoundary)};
    nodesBoundary = getNodes(data.mesh.(nameConBoundary));
    data.ip.waveNumberValue = data.bottom.waveNumber(nodesBoundary(1));
end
angle = data.ip.direction*(pi/180);
kvector = data.ip.waveNumberValue*[cos(angle) sin(angle)];
[data.ip.value,data.ip.gradientsIP,data.mesh.der2DShapeFunOn1D,data.mesh.extT,...
    data.mesh.extTb,data.mesh.intT,data.mesh.inPMLelems,data.mesh.elementFaceInfo.boundaryPML] =...
    computeIncidentWaveField(data,kvector,omega,handle);

% Global absorption PML parameter and stretching
data.PMLabsorptionValue = zeros(size(data.mesh.(nameNodal)));
data.PMLstretching = data.mesh.(nameNodal);
for iboundary = 1:numel(data.mesh.boundaryNames)
    if strcmp(data.PML{5,iboundary},'on')
        data.PMLabsorptionValue = data.PMLabsorptionValue + ...
            data.PML{1,iboundary};
        %data.PMLabsorptionValue = 0 * data.PMLabsorptionValue; %%%%%%%%%%%%%%%%%%%%
        data.PMLstretching = data.PMLstretching + data.PML{7,iboundary}; 
    end
end

% Volume matrices and volume vector from Berkhoff equation
setOutput({'Computing stiffness, mass matrices and volume vector (FEM)...'},handle)
t = cputime;
tic
if data.staticCondensation
    setOutput({'Static condensation enabled'},handle)
    [data.mesh.MminusK,data.mesh.fvolume,data.mesh.elementalInnerInfo,data.mesh.newFaceNodesNum] =...
        berkhoffMatricesAndVolumeVector_static(data.mesh.(nameNodal),...
        data.mesh.(nameCon),...
        data.mesh.referenceElement,...
        data.ip.amplitude,...
        kvector,...
        data.bottom.waveNumber,...
        data.bottom.ccg,...
        data.PMLabsorptionValue,...
        omega);
else
    [data.mesh.MminusK,data.mesh.fvolume] = berkhoffMatricesAndVolumeVector...
        (data.mesh.(nameNodal),...
        data.mesh.(nameCon),...
        data.mesh.referenceElement,...
        data.ip.value,...
        data.bottom.waveNumber,...
        data.bottom.ccg,...
        data.PMLabsorptionValue,...
        omega,...
        data.mesh.inPMLelems);
end

data.cputime.volMatToc = toc;
data.cputime.volMat = cputime - t;



