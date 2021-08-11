function BC_PMLvalue_R_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = str2double(get(hObject,'String'));
if isnan(value) || value < 0
    set(handles.BC_PMLvalue_R,'String','value')
    set(handles.BC_set,'Enable','off')
    msje = {'??? The R parameter has to be a positive or zero number'};
    setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
end
data.currentBoundary.parameter(2) = value;
param = data.currentBoundary.parameter;
set(handles.BC_set,'Enable','on')
for cont = [1 3]
    if param(cont) == -1
        set(handles.BC_set,'Enable','off')
        break
    end
end
       
guidata(hObject,data);