function createMeshData(file,ext,handles)

data = guidata(handles.MainFigure);

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
    auxMeshFields = fieldnames(auxMesh);
    boundary = {};
    submesh = {};
    indexBoundary = [];
    indexSubmesh = [];
    boundary = {};
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
    nOfNodes = size(auxMesh.(auxMeshFields{indexElemPosCon(2)}),1);
    nOfElements = size(auxMesh.(auxMeshFields{indexElemPosCon(3)}),1);

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
    nOfNodes = size(auxMesh.X,1);
    nOfElements = size(auxMesh.T,1) + sizeBoundary;
    k = length(indexBoundary) + 1;
    u = length(indexSubmesh) + 1;
end

clear auxMesh

set(handles.BC_PMLvalue_area,'String',{'none' submesh{:}},'Value',1)

setOutput({['Detected ' num2str(k-1) ' boundaries and ' num2str(u-1) ' subdomains']
    ['Nodes: ' num2str(nOfNodes) '   Elements: ' num2str(nOfElements)]},handles.run_wipOutput)

% Condition for plotting
data.plotCondition = nOfNodes < 1.5e6 && nOfElements < 1e6;
if ~data.plotCondition
    setOutput({'Mesh size overpassed'
               'Plots are not going to be visualized'},handles.run_wipOutput)
end

data.staticCondensation = false;
data.constantBottomFlag = false;

guidata(handles.MainFigure,data);


%_________________
    function findBoundaryAndCheckFEM
        
        elemInfoCheck = false;
        nodalPosCheck = false;
        conectivityCheck = false;
        typeCheck = false;
        nOfNodesCheck = false;
        faces2dCheck = false;
        faces1dCheck = false;
        elementFaceInfoCheck = false;
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
            elseif strcmpi(iboundary,'elementFaceInfo')
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

        if warnFlag
            msje = {['??? The mesh might has been well read, but some others unnecessary '...
                'variables have been detected too. To avoid problems please load a '...
                'mesh file which contains only necessary variables. Check the '...
                'Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        elseif ~numel(boundary)
            msje = {['??? The boundary cant be read. Make sure that the file contains '...
                'the conectivity matrices "T_boundaryName" for the boundary '...
                '"boundaryName". Check the Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        elseif ~elemInfoCheck
            msje = {['??? The element info cant be read. Make sure that the file '...
                'contains a structure which has the correct element info. '...
                'Check the Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        elseif ~typeCheck
            msje = {['??? The element type cant be read. Make sure that the element '...
                'info contains a scalar 0 for QUA type or 1 for TRI type. Check the '...
                'Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        elseif ~nOfNodesCheck
            msje = {['??? The referent element�s number of nodes cant be read. '...
                'Make sure that the element info contains a scalar with the number '...
                'of nodes. Check the Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        elseif ~faces2dCheck
            msje = {['??? The 2D face�s element nodes info cant be read. Make sure that '...
                'the element info contains a matrix with the nodes� numbering for '...
                'each face. Check the Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        elseif ~faces1dCheck
            msje = {['??? The 1D element nodes info cant be read. Make sure that '...
                'the element info contains a vector with the 1D nodes� numbering. '...
                'Check the Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        elseif ~nodalPosCheck
            msje = {['??? The nodal position matrix cant be read. Make sure that '...
                'the file contains a correct nodal position matrix. Check '...
                'the Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        elseif ~conectivityCheck
            msje = {['??? The conectivity matrix cant be read. Make sure that '...
                'the file contains a correct conectivity matrix. Check '...
                'the Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        elseif strcmp(data.computation,'CDG') && ~elementFaceInfoCheck
            msje = {['??? The face element info cant be read. Make sure that the file '...
                'contains a structure which has the correct face element info. '...
                'Check the Berkhoff GUI help for details']};
            setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
        end

        indexElemPosCon = [indexElemInfo indexPosMatrix indexConectivity];
        indexTypNodFac1Fac2 = [indexType nOfNodesIndex faces1dIndex faces2dIndex];

    end
%_________________


end
