function run_all_Callback(hObject,eventdata,handles)

if get(hObject,'Value')
    set([handles.run_domain handles.run_boundary],'Value',0)
    set([handles.res_display handles.res_save],'Enable','off')

    data = guidata(handles.MainFigure);


    setOutput({'-----------------------------------------'
        'Work in process...'},handles.run_wipOutput)
    if strcmp(data.computation,'FEM')
        data = run_domain(data,handles.run_wipOutput);
        data = run_boundary(data,handles.run_wipOutput);
%         data = run_domain_REFLECTEDWAVE_OLD(data,handles.run_wipOutput);
%         data = run_boundary_REFLECTEDWAVE_OLD(data,handles.run_wipOutput);
    elseif strcmp(data.computation,'NEFEM')
        data = run_domain_NEFEM(data,handles.run_wipOutput);
        data = run_boundary_NEFEM(data,handles.run_wipOutput);
    elseif strcmp(data.computation,'CDG')
        data = run_domain_CDG(data,handles.run_wipOutput);
        data = run_boundary_CDG(data,handles.run_wipOutput);
    elseif strcmp(data.computation,'DG')
        data = preprocess_DG(data,handles.run_wipOutput);
        data = run_simulation_DG(data,handles.run_wipOutput);
    end
    set([handles.res_display handles.res_save handles.run_boundary],...
        'Enable','on')
    setOutput({'done!'},handles.run_wipOutput)

    guidata(handles.MainFigure,data);
else
    set([handles.res_display handles.run_boundary],'Enable','off')
end
