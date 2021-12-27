function ip_amplitudeValue_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = str2double(get(hObject,'String'));
if isnan(value) || value <= 0
    if ~isempty(data.ip.amplitude)
        set(handles.ip_amplitudeValue,'String',num2str(data.ip.amplitude))
    else
        set(handles.ip_amplitudeValue,'String','value')
    end
    msje = {'??? The amplitude has to be a positive number'};
    setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
elseif isempty(data.ip.amplitude) || value ~= data.ip.amplitude
    data.ip.amplitude = value;
    if isfield(data.ip,'value')
        data.ip = rmfield(data.ip,'value');
    end
    
    guidata(hObject,data);
    
    check_run(handles)
end