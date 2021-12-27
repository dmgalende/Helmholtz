function panXY_Callback(hObject,eventdata,handles)

switch get(hObject,'State')
    case 'on'
        set(handles.zoomXY,'State','off')
        pan(handles.axesHandle,'on')
    case 'off'
        pan(handles.axesHandle,'off')
end