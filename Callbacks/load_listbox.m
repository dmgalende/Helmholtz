function load_listbox(dir_path,handle,fieldName)

data = guidata(handle(1));

for i = 1:numel(dir_path)
    ipath = dir_path{i};
    dir_aux = dir(ipath);
    [sorted_names,sorted_index] = sortrows({dir_aux.name}');
    data.dir_struct.(fieldName{i}).files_names = sorted_names;
    data.dir_struct.(fieldName{i}).sorted_index = sorted_index;
    data.dir_struct.(fieldName{i}).isdir = [dir_aux.isdir];
    data.dir_struct.(fieldName{i}).cd = ipath;
    set(handle(i),'String',sorted_names,'Value',1)
end

guidata(handle(1),data);
    
    
    
