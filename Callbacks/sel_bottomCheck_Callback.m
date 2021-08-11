function sel_bottomCheck_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = get(hObject,'Value');
if value
    set(data.plotBottom,'Visible','on')
else
    set(data.plotBottom,'Visible','off')
end

guidata(hObject,data);