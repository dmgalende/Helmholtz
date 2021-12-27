function mesh_list_Callback(hObject,eventdata,handles)

if strcmp(get(handles.MainFigure,'SelectionType'),'open')
    data = guidata(hObject);
    dir_struct = data.dir_struct;
    
    index_selected = get(hObject,'Value');
    file_list = get(hObject,'String');
    filename = file_list{index_selected};
    directory = dir_struct.mesh.cd;
    if dir_struct.mesh.isdir(dir_struct.mesh.sorted_index(index_selected))
        dir_path = [directory '/' filename];
        load_listbox({dir_path},hObject,{'mesh'})
    else
        [path,name,ext] = fileparts(filename);
        cd(directory)
        
        setOutput({'Loading mesh data...'},handles.run_wipOutput,1)
        createMeshData(filename,ext,handles)
        
        % Computing reference element and boundary, mesh plots
        elemType = getReferenceElement(handles);
        setOutput({'Initializing plots...'},handles.run_wipOutput)
        if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
            initializePlots(handles)
        elseif strcmp(data.computation,'NEFEM')
            initializePlots_NEFEM(handles);
        end

        resetBC(handles)
        set(handles.sel_mesh,'String',filename)
        set(handles.sel_elem,'String',elemType)
        set(handles.sel_bottom,'String','none')
        handles2Enable = [
            handles.bottom_list 
            handles.sel_meshCheck
            handles.ip_boundaryCheck
            handles.BC_list
            handles.BC_popup
            handles.BC_set
            handles.zoomXY
            handles.panXY
            handles.compOptions
         	];
        handles2Disable = [
            handles.run_domain
            handles.run_boundary
            handles.run_all
            handles.sel_bottomCheck
            handles.res_display
            handles.res_save
            ];
        set(handles2Enable,'Enable','on')
        set(handles2Disable,'Enable','off','Value',0)
        set(handles.fixBottom,'Enable','off')
        if strcmp(data.computation,'DG')
            set(handles.SimParam,'Enable','on')
        else
            set(handles.SimParam,'Enable','off')
        end
        setOutput({'done!'},handles.run_wipOutput)
    end
end
