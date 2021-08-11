function run_boundary_Callback(hObject,eventdata,handles)

if get(hObject,'Value')

    data = guidata(handles.MainFigure);


    setOutput({'-----------------------------------------'
        'Work in process...'},handles.run_wipOutput)
    if strcmp(data.computation,'FEM')
        data = run_boundary(data,handles.run_wipOutput);
    elseif strcmp(data.computation,'NEFEM')
        data = run_boundary_NEFEM(data,handles.run_wipOutput);
    elseif strcmp(data.computation,'CDG')
        data = run_boundary_CDG(data,handles.run_wipOutput);
    elseif strcmp(data.computation,'DG')
        data = run_simulation_DG(data,handles.run_wipOutput);
    end
    set([handles.res_display handles.res_save],'Enable','on')
    setOutput({'done!'},handles.run_wipOutput)

    guidata(handles.MainFigure,data);
else
    set(handles.res_display,'Enable','off')
end
