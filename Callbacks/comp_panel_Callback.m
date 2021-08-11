function comp_panel_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'comp_FEM'
        data.computation = 'FEM';
        set(handles.run_domain,'String','Run domain')
        set(handles.run_boundary,'String','Run boundary')
    case 'comp_NEFEM'
        data.computation = 'NEFEM';
        set(handles.run_domain,'String','Run domain')
        set(handles.run_boundary,'String','Run boundary')
    case 'comp_CDG'
        data.computation = 'CDG';
        set(handles.run_domain,'String','Run domain')
        set(handles.run_boundary,'String','Run boundary')
    case 'comp_DG'
        data.computation = 'DG';
        set(handles.run_domain,'String','Preprocess')
        set(handles.run_boundary,'String','Run')
end

guidata(handles.MainFigure,data); 
