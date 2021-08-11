function res_load_Callback(hObject,eventdata,handles)

data = guidata(hObject);

[fileName,path] = uigetfile('*.mat','Load data');
if ischar(fileName) && ~strcmp(fileName(end:-1:end-3),'tam.')
    setOutput({['??? Error loading data file. The file extension has to '...
        'be .mat']},handles.run_wipOutput), error('Error in Berkhoff GUI')
elseif fileName == 0
    setOutput({'Load canceled'},handles.run_wipOutput)
    return
end

% Variable to preserve
dir_struct = data.dir_struct;

% Load data
setOutput({'Loading data file...'},handles.run_wipOutput,1)
clear data
loadVar = load([path '/' fileName]);
data = loadVar.data;
data.dir_struct = dir_struct;

% For older GUI versions
if ~isfield(data,'plotCondition'), data.plotCondition = true; end
if ~isfield(data,'computation'), data.computation = 'FEM'; end
if ~isfield(data,'staticCondensation'), data.staticCondensation = false; end
if ~isfield(data,'constantBottomFlag'), data.constantBottomFlag = false; end
if size(data.PML,1) < 7, data.PML(7,:) = {0}; end

guidata(hObject,data);

% INITIALIZATION TASKS

% Plots and selection info
if ~data.plotCondition
    setOutput({'Mesh size overpassed'
        'Plots are not going to be visualized'},handles.run_wipOutput)
end
setOutput({'Initializing plots...'},handles.run_wipOutput)
if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
    nameElem = data.mesh.fieldNames{data.mesh.indexElemPosCon(1)};
    nameType = data.mesh.elemFieldNames{data.mesh.indexTypNodFac1Fac2(1)};
    nameNodes = data.mesh.elemFieldNames{data.mesh.indexTypNodFac1Fac2(2)};
    if data.mesh.(nameElem).(nameType)
        elemString = 'Tri ';
    else
        elemString = 'Qua ';
    end
    elemType = [elemString 'with ' num2str(data.mesh.(nameElem).(nameNodes))...
        ' nodes'];
    if strcmp(data.computation,'FEM')
        set(handles.comp_FEM,'value',1)
    elseif strcmp(data.computation,'CDG')
        set(handles.comp_CDG,'value',1)
    else
        set(handles.comp_DG,'value',1)
    end
    initializePlots(handles)
elseif strcmp(data.computation,'NEFEM')
    elemType = ['Tri with ' num2str(data.mesh.elemInfo.nOfNodes) ' nodes'];
    set(handles.comp_NEFEM,'value',1)
    initializePlots_NEFEM(handles)
end
plotBottom(handles)
set([handles.sel_mesh handles.sel_bottom],'String',['Loaded from ' fileName])
set(handles.sel_elem,'String',elemType)

data = guidata(hObject);

% Enable all commands
handlesNames = fieldnames(handles);
for iset = 1:numel(handlesNames)
    ihandle = handles.(handlesNames{iset});
    typeHandle = get(ihandle,'Type');
    if ~any(strcmp(typeHandle,{'figure' 'uipanel' 'axes' 'uibuttongroup'}))
        set(ihandle,'Enable','on')
    end
end

% Incident potential data
setOutput({'Recovering incident potential data...'},handles.run_wipOutput)
set(handles.ip_waveNumberBoundary,'String',{'none' data.mesh.boundaryNames{:}},...
    'Value',data.ip.waveNumberBoundary+1)
if data.ip.waveNumberBoundary
    set(handles.ip_boundaryCheck,'Value',1)
    set(handles.ip_waveNumberBoundary,'Visible','on')
    set(handles.ip_waveNumberValue,'Visible','off')
else
    set(handles.ip_boundaryCheck,'Value',0)
    set(handles.ip_waveNumberBoundary,'Visible','off')
    set(handles.ip_waveNumberValue,'Visible','on')
    set(handles.ip_waveNumberValue,'String',num2str(data.ip.waveNumberValue))
end
set(handles.ip_directionValue,'String',num2str(data.ip.direction))
set(handles.ip_amplitudeValue,'String',num2str(data.ip.amplitude))
set(handles.ip_periodValue,'String',num2str(data.ip.period))

% Boundary conditions
setOutput({'Recovering BC data...'},handles.run_wipOutput)
set(handles.BC_list,'String',data.mesh.boundaryNames,'Value',1)
set(handles.BC_PMLvalue_area,'String',{'none' data.mesh.submeshNames{:}})
data.currentBoundary.value = 1;
data.currentBoundary.condition = data.BC.values(1);
set(handles.BC_popup,'Value',data.currentBoundary.condition)

% Boundary Conditions' plots
setOutput({'Initializing BC plots...'},handles.run_wipOutput)
nOfBoundaries = numel(data.mesh.boundaryNames);
data.plotPMLarea(1,:) = cell(1,nOfBoundaries);
data.plotPMLnodes = cell(2,nOfBoundaries);
if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
    nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
elseif strcmp(data.computation,'NEFEM')
    nameNodal = 'X';
end
for aux = 1:nOfBoundaries
    if data.BC.values(aux) == 4
        param = data.BC.parameters{aux}{4};

        % PML area plot
        data.plotPMLarea{2,aux} = param(3);
        data.plotPMLarea{1,aux} = plotPMLarea(data.PML{2,aux},param(3),handles.axesHandle);

        % PML nodes plot
        if ~isempty(data.PML{6,aux})
            PMLnodes = data.PML{6,aux};
        else
            PMLnodes = 1:size(data.mesh.(nameNodal),1);
        end
        data.plotPMLnodes(:,aux) = ...
            plotPMLnodes(data.mesh.(nameNodal)(PMLnodes,:),data.PML{4,aux},handles.axesHandle);
    end
end

%Checking if the solution has been computed
if isfield(data,'solution')
    if any(strcmp(data.computation,{'FEM' 'NEFEM'}))
        cond = all(isfield(data.mesh,{'MminusK' 'fvolume'}));
    elseif strcmp(data.computation,'CDG')
        cond = all(isfield(data.mesh,{'KminusMplusH' 'fvolume'}));
    end
    if ~cond, set(handles.run_boundary,'Enable','off'), end
else
    set([handles.res_display handles.run_boundary],'Enable','off')
end

guidata(hObject,data);

check_BC(data.currentBoundary.condition,handles)

set([handles.run_all handles.run_boundary handles.run_domain],'Value',0)
setOutput({'done!'},handles.run_wipOutput)















