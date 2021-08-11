function BC_PMLareaCheck_Callback(hObject,eventdata,handles)

data = guidata(hObject);

if get(hObject,'Value')
    set(data.plotPMLarea{1,data.currentBoundary.value},'Visible','on')
else
    set(data.plotPMLarea{1,data.currentBoundary.value},'Visible','off')
end