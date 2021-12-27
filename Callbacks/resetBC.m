function resetBC(handles)

data = guidata(handles.MainFigure);
if isfield(data,'BC')
    data = rmfield(data,{'BC' 'currentBoundary'});
end

set(handles.BC_list,'String',data.mesh.boundaryNames)
set(handles.ip_waveNumberBoundary,'String',{'none' data.mesh.boundaryNames{:}},...
    'Value',1);
data.ip.waveNumberBoundary = 0;
data.BC = struct;
nOfBoundaries = numel(data.mesh.boundaryNames);
nOfConditions = numel(get(handles.BC_popup,'String'));
data.BC.values = ones(nOfBoundaries,1);
data.BC.parameters = cell(nOfBoundaries,1);
data.BC.parameters(:) = {cell(nOfConditions,1)};
data.PML = cell(7,nOfBoundaries);
data.PML(5,:) = {'off'};
data.plotPMLarea = cell(2,nOfBoundaries);
data.plotPMLnodes = cell(2,nOfBoundaries);
data.currentBoundary = struct('value',1,'condition',1,'parameter',[]);
set(handles.BC_list,'Value',data.currentBoundary.value)
set(handles.BC_popup,'Value',data.currentBoundary.condition)

guidata(handles.MainFigure,data);

check_BC(data.currentBoundary.value,handles)