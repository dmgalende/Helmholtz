function bottom_list_Callback(hObject,eventdata,handles)

if strcmp(get(handles.MainFigure,'SelectionType'),'open')
    data = guidata(hObject);
    dir_struct = data.dir_struct;
    
    index_selected = get(hObject,'Value');
    file_list = get(hObject,'String');
    filename = file_list{index_selected};
    directory = dir_struct.bottom.cd;
    if dir_struct.bottom.isdir(dir_struct.bottom.sorted_index(index_selected))
        dir_path = [directory '/' filename];
        load_listbox({dir_path},hObject,{'bottom'})
    else
        [path,name,ext] = fileparts(filename);
        cd(directory)
        
        setOutput({'Loading bottom data...'},handles.run_wipOutput)
        flag = createBottomData(filename,ext,handles);
        if ~flag
            setOutput({'Canceled by user'},handles.run_wipOutput)
            return
        end
        plotBottom(handles)
        set(handles.sel_bottom,'String',filename)
        set([handles.sel_bottomCheck handles.fixBottom],'Enable','on')
        set([handles.run_domain handles.run_all],'Value',0)
        check_run(handles)
        setOutput({'done!'},handles.run_wipOutput)
    end
end