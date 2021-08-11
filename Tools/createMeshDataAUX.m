function data = createMeshDataAUX(data,file)

ext = file(end-3:end);

if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
    switch ext
        case '.mat'
            auxMesh = load(file);
        case '.dcm'
            if strcmp(data.computation,'FEM')
                meshfile = GenerateMatFileFromEZ4U(file);
            elseif strcmp(data.computation,'CDG')
                meshfile = GenerateMatFileFromEZ4U_CDG(file);
            else
                meshfile = GenerateMatFileFromEZ4U_DG(file);
            end
            auxMesh = load(meshfile{:});
            delete(meshfile{:})
        otherwise
            msje = {['??? The mesh has to be loaded from a .mat file or a .dcm EZ4U '...
                'mesh file. Check the Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
    end
    findBoundaryAndCheckFEM

    if isfield(data,'mesh')
        data = rmfield(data,'mesh');
        if isfield(data,'bottom')
            data = rmfield(data,'bottom');
        end
        if isfield(data.ip,'value')
            data.ip = rmfield(data.ip,'value');
        end
        if isfield(data,'infoFaces')
            data = rmfield(data,'infoFaces');
        end
    end
    if isfield(data,'solution')
        data = rmfield(data,'solution');
    end
    data.mesh = auxMesh;
    data.mesh.fieldNames = auxMeshFields;
    data.mesh.boundaryNames = boundary;
    data.mesh.boundaryIndex = indexBoundary;
    data.mesh.submeshNames = submesh;
    data.mesh.submeshIndex = indexSubmesh;
    data.mesh.indexElemPosCon = indexElemPosCon;
    data.mesh.elemFieldNames = subfields;
    data.mesh.indexTypNodFac1Fac2 = indexTypNodFac1Fac2;

elseif strcmp(data.computation,'NEFEM')
    switch ext
        case '.mat'
            auxMesh = load(file);
        otherwise
            msje = {['??? For NEFEM computation, the mesh has to be loaded '...
                'from the .mat file done by NEFEM preprocess']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
    end

    if isfield(data,'mesh')
        data = rmfield(data,'mesh');
        if isfield(data,'bottom')
            data = rmfield(data,'bottom');
        end
        if isfield(data.ip,'value')
            data.ip = rmfield(data.ip,'value');
        end
    end
    if isfield(data,'solution')
        data = rmfield(data,'solution');
    end
    auxMeshFields = fieldnames(auxMesh);
    boundary = {};
    submesh = {};
    indexBoundary = [];
    indexSubmesh = [];
    sizeBoundary = 0;
    for i = 1:length(auxMeshFields)
        name = auxMeshFields{i};
        if length(name) > 3 && strcmpi(name(1:3),'Tb_')
            indexBoundary = [indexBoundary i];
            sizeBoundary = sizeBoundary + size(auxMesh.(name),1);
            boundary = {boundary{:} name(4:end)};
        elseif length(name) > 2 && strcmpi(name(1:2),'T_')
            indexSubmesh = [indexSubmesh i];
            submesh = {submesh{:} name(3:end)};
        end
    end
    data.mesh = auxMesh;
    data.mesh.fieldNames = auxMeshFields;
    data.mesh.boundaryNames = boundary;
    data.mesh.boundaryIndex = indexBoundary;
    data.mesh.submeshNames = submesh;
    data.mesh.submeshIndex = indexSubmesh;
end


%_________________
    function findBoundaryAndCheckFEM

        auxMeshFields = fieldnames(auxMesh);
        boundary = {};
        submesh = {};
        indexBoundary = [];
        indexSubmesh = [];
        elemInfoCheck = false;
        nodalPosCheck = false;
        conectivityCheck = false;
        typeCheck = false;
        nOfNodesCheck = false;
        faces2dCheck = false;
        faces1dCheck = false;
        elementFaceInfoCheck = false; %Only for CDG
        warnFlag = false;
        k = 1; %index for boundary
        u = 1; %index for submeshes

        for ifield = 1:numel(auxMeshFields)
            iboundary = auxMeshFields{ifield};
            idata = auxMesh.(iboundary);
            if length(iboundary) > 3 && strcmpi(iboundary(1:3),'Tb_') &&...
                    isnumeric(idata) && all(size(idata) >= [1 2]) %Optimized for CDG
                boundary{k} = iboundary(4:end);
                indexBoundary = [indexBoundary ifield];
                k = k + 1;
            elseif strcmpi(iboundary,'elementFaceInfo') %Only for CDG
                elementFaceInfoCheck = true;
            elseif ~elemInfoCheck && isstruct(idata)
                subfields = fieldnames(idata);
                for strField = 1:numel(subfields)
                    subdata = idata.(subfields{strField});
                    if ~typeCheck && isnumeric(subdata) && isscalar(subdata)...
                            && any(subdata == [0 1])
                        typeCheck = true;
                        indexType = strField;
                    elseif ~nOfNodesCheck && isnumeric(subdata) &&...
                            isscalar(subdata) && subdata >= 3
                        nOfNodesCheck = true;
                        nOfNodesIndex = strField;
                    elseif ~faces2dCheck && isnumeric(subdata) &&...
                            all(size(subdata) >= [3 2])
                        faces2dCheck = true;
                        faces2dIndex = strField;
                    elseif ~faces1dCheck && isnumeric(subdata) &&...
                            any(size(subdata) == 1) && length(subdata) >= 2
                        faces1dCheck = true;
                        faces1dIndex = strField;
                    elseif ~elementFaceInfoCheck
                        warnFlag = true;
                        break
                    end
                end
                elemInfoCheck = true;
                indexElemInfo = ifield;
                if warnFlag, break, end
            elseif ~nodalPosCheck && isnumeric(idata) &&...
                    size(idata,1) >= 3 && size(idata,2) == 2
                nodalPosCheck = true;
                indexPosMatrix = ifield;
            elseif isnumeric(idata) && all(size(idata) >= [1 3])
                if ~conectivityCheck
                    conectivityCheck = true;
                    indexConectivity = ifield;
                elseif length(iboundary) > 2 && strcmpi(iboundary(1:2),'T_')
                    submesh{u} = iboundary(3:end);
                    indexSubmesh = [indexSubmesh ifield];
                    u = u + 1;
                else
                    warnFlag = true;
                    break
                end
            end
        end

        indexElemPosCon = [indexElemInfo indexPosMatrix indexConectivity];
        indexTypNodFac1Fac2 = [indexType nOfNodesIndex faces1dIndex faces2dIndex];

    end
%_________________


end