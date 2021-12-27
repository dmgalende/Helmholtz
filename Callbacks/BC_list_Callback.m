function BC_list_Callback(hObject,eventdata,handles)

data = guidata(hObject);
    
set(data.plotBoundary{data.currentBoundary.value},'Color','r')
data.currentBoundary.value = get(hObject,'Value');

guidata(hObject,data);

set(data.plotBoundary{data.currentBoundary.value},'Color','g')
set(handles.BC_popup,'Value',data.BC.values(data.currentBoundary.value))
check_BC(data.BC.values(data.currentBoundary.value),handles)

