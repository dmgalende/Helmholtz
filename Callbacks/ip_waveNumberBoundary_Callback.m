function ip_waveNumberBoundary_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = get(hObject,'Value') - 1;
if value ~= data.ip.waveNumberBoundary
    data.ip.waveNumberBoundary = value;
    if isfield(data.ip,'value')
        data.ip = rmfield(data.ip,'value');
    end
    
    guidata(hObject,data);
    
    check_run(handles)
end
