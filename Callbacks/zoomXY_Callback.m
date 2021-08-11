function zoomXY_Callback(hObject,eventdata,handles)

switch get(hObject,'State')
    case 'on'
        set(handles.panXY,'State','off')
        zoom(handles.axesHandle,'on')
        zoomHandle = zoom(handles.axesHandle);
        set(zoomHandle,'RightClickAction','InverseZoom')
    case 'off'
        zoom(handles.axesHandle,'off')
end
        