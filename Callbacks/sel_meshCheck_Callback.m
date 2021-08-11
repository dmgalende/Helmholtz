function sel_meshCheck_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = get(hObject,'Value');
if value
    set(data.plotMesh,'Visible','on')
else
    set(data.plotMesh,'Visible','off')
end

guidata(hObject,data);