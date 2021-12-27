function ip_directionValue_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = str2double(get(hObject,'String'));
if isnan(value)
    if ~isempty(data.ip.direction)
        set(handles.ip_directionValue,'String',num2str(data.ip.direction))
    else
        set(handles.ip_directionValue,'String','value')
    end
    msje = {'??? The direction has to be a number'};
    setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
elseif isempty(data.ip.direction) || value ~= data.ip.direction
    data.ip.direction = value;
    set([handles.run_domain handles.run_all],'Value',0)
    set([handles.run_boundary handles.res_save handles.res_display],...
        'Enable','off')
    if isfield(data.ip,'value')
        data.ip = rmfield(data.ip,'value');
    end
    
    guidata(hObject,data);
    
    check_run(handles)
end

