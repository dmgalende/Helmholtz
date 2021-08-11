function BC_PMLvalue_area_Callback(hObject,eventdata,handles)

data = guidata(hObject);

value = get(hObject,'Value');
data.currentBoundary.parameter(4) = value;

if numel(data.mesh.submeshNames) > 0
    set(data.plotSubmesh,'Visible','off')
    if value > 1
        set(data.plotSubmesh(value-1),'Visible','on')
    end
end

guidata(hObject,data);