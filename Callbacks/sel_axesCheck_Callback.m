function sel_axesCheck_Callback(hObject,eventdata,handles)

value = get(hObject,'Value');
if value
    set(handles.axesHandle,'Visible','on')
else
    set(handles.axesHandle,'Visible','off')
end