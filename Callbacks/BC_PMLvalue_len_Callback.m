function BC_PMLvalue_len_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = str2double(get(hObject,'String'));
if isnan(value) || value <= 0 
    set(handles.BC_PMLvalue_len,'String','value')
    set(handles.BC_set,'Enable','off')
    msje = {'??? The length of the PML region has to be a positive number'};
    setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
end
data.currentBoundary.parameter(3) = value;
param = data.currentBoundary.parameter;
set(handles.BC_set,'Enable','on')
for cont = [1 2]
    if param(cont) == -1
        set(handles.BC_set,'Enable','off')
        break
    end
end

%Plot user selected PML area
data.plotPMLarea{2,data.currentBoundary.value} = value;
delete(data.plotPMLarea{1,data.currentBoundary.value})          
data.plotPMLarea{1,data.currentBoundary.value} = ...
    plotPMLarea(data.PML{2,data.currentBoundary.value},value,handles.axesHandle);

%Visibility
if ~get(handles.BC_PMLareaCheck,'Value')
    set(data.plotPMLarea{1,data.currentBoundary.value},'Visible','off')
end

set(handles.BC_PMLareaCheck,'Enable','on')
                            
guidata(hObject,data);
