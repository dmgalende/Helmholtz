function ip_periodValue_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = str2double(get(hObject,'String'));
if isnan(value) || value <= 0
    if ~isempty(data.ip.period)
        set(handles.ip_periodValue,'String',num2str(data.ip.period))
    else
        set(handles.ip_periodValue,'String','value')
    end
    msje = {'??? The period has to be a positive number'};
    setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
elseif isempty(data.ip.period) || value ~= data.ip.period
    data.ip.period = value;
    if isfield(data,'bottom') && isfield(data.bottom,'ccg')
        data.bottom = rmfield(data.bottom,{'waveNumber' 'ccg'});
    end
    if isfield(data.ip,'value')
        data.ip = rmfield(data.ip,'value');
    end
    
    guidata(hObject,data);
    
    check_run(handles)    
end

