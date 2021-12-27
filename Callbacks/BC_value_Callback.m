function BC_value_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = str2double(get(hObject,'String'));
if data.currentBoundary.condition == 2 && (isnan(value) || value < 0 ||...
        value > 1)
    set(handles.BC_value,'String','value')
    set(handles.BC_set,'Enable','off')
    msje = {['??? The alpha parameter has to be a number between 0 and 1, both '...
        'values included']};
    setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
elseif data.currentBoundary.condition == 3 && (isnan(value) || value <= 0)
    set(handles.BC_value,'String','value')
    set(handles.BC_set,'Enable','off')
    msje = {'??? The radius has to be a positive number'};
    setOutput(msje,handles.run_wipOutput), error('Error in Berkhoff GUI')
end
data.currentBoundary.parameter = value;
set(handles.BC_set,'Enable','on')

guidata(hObject,data);