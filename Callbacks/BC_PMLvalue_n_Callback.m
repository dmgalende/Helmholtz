function BC_PMLvalue_n_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = str2double(get(hObject,'String'));
if isnan(value) || value < 0 
    set(handles.BC_PMLvalue_n,'String','value')
    set(handles.BC_set,'Enable','off')
    msje = {'??? The absorption degree has to be a positive or zero number'};
    setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
end
data.currentBoundary.parameter(1) = value;
param = data.currentBoundary.parameter;
set(handles.BC_set,'Enable','on')
for cont = [2 3]
    if param(cont) == -1
        set(handles.BC_set,'Enable','off')
        break
    end
end
       
guidata(hObject,data);