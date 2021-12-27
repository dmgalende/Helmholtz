function ip_waveNumberValue_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = str2double(get(hObject,'String'));
if isnan(value) || value <= 0
    if ~isempty(data.ip.waveNumberValue)
        set(handles.ip_waveNumberValue,'String',num2str(data.ip.waveNumberValue))
    else
        set(handles.ip_waveNumberValue,'String','value')
    end
    msje = {'??? The wave number has to be a positive number'};
    setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
elseif isempty(data.ip.waveNumberValue) || value ~= data.ip.waveNumberValue
    data.ip.waveNumberValue = value;
    if isfield(data.ip,'value')
        data.ip = rmfield(data.ip,'value');
    end
    
    guidata(hObject,data);
    
    check_run(handles)
end

