function ip_boundaryCheck_Callback(hObject,eventdata,handles)

value = get(hObject,'Value');
if value
    set(handles.ip_waveNumberBoundary,'Visible','on')
    set(handles.ip_waveNumberValue,'Visible','off')
else
    data = guidata(hObject);
    
    set(handles.ip_waveNumberBoundary,'Visible','off','Value',1)
    data.ip.waveNumberBoundary = 0;
    data.ip.waveNumberValue = [];
    set(handles.ip_waveNumberValue,'Visible','on','String','value')
    
    guidata(hObject,data);
end

check_run(handles)