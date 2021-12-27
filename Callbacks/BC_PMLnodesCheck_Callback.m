function BC_PMLnodesCheck_Callback(hObject,eventdata,handles)

data = guidata(hObject);

if get(hObject,'Value')
    set([data.plotPMLnodes{:,data.currentBoundary.value}],'Visible','on')
else
    set([data.plotPMLnodes{:,data.currentBoundary.value}],'Visible','off')
end