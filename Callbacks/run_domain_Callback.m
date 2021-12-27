function run_domain_Callback(hObject,eventdata,handles)

if get(hObject,'Value')
    set([handles.run_all handles.run_boundary],'Value',0)
    set([handles.res_display handles.res_save],'Enable','off')

    data = guidata(handles.MainFigure);


    setOutput({'-----------------------------------------'
        'Work in process...'},handles.run_wipOutput)
    if strcmp(data.computation,'FEM')
        data = run_domain(data,handles.run_wipOutput);
    elseif strcmp(data.computation,'NEFEM')
        data = run_domain_NEFEM(data,handles.run_wipOutput);
    elseif strcmp(data.computation,'CDG')
        data = run_domain_CDG(data,handles.run_wipOutput);
    elseif strcmp(data.computation,'DG')
        data = preprocess_DG(data,handles.run_wipOutput);
        set(handles.res_save,'Enable','on')
    end
    set(handles.run_boundary,'Enable','on')
    setOutput({'done!'},handles.run_wipOutput)

    guidata(handles.MainFigure,data);
else
    set(handles.run_boundary,'Enable','off','Value',0)
end

