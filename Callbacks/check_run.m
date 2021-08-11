function check_run(handles)

data = guidata(handles.MainFigure);

period = data.ip.period;
direction = data.ip.direction;
amplitude = data.ip.amplitude;
waveNumberCondition = getWaveNumberCondition;
handles2Change = [handles.run_domain handles.run_all handles.res_save];
if isfield(data,'bottom') && ~isempty(period) && ~isempty(direction)...
       && ~isempty(amplitude) && waveNumberCondition
    set(handles2Change,'Enable','on','Value',0)
else
    set(handles2Change,'Enable','off','Value',0)
end
set([handles.run_boundary handles.res_display],'Enable','off','Value',0)
if isfield(data,'solution')
    data = rmfield(data,'solution');
    
    guidata(handles.MainFigure,data);
end

    %_________________
    function isIt = getWaveNumberCondition
    
    waveNumberFromBoundary = get(handles.ip_boundaryCheck,'Value');
    if waveNumberFromBoundary && data.ip.waveNumberBoundary
        isIt = true;
    elseif ~waveNumberFromBoundary && ~isempty(data.ip.waveNumberValue)
        isIt = true;
    else
        isIt = false;
    end
    
    end
    %_________________


end