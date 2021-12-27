function res_save_Callback(hObject,eventdata,handles)

[fileName,path] = uiputfile('*.mat','Save data as');
if ischar(fileName) && ~strcmp(fileName(end:-1:end-3),'tam.')
    setOutput({['??? Error saving the file. The file extension has to '...
        'be .mat']},handles.run_wipOutput), error('Error in Berkhoff GUI')
elseif fileName == 0
    setOutput({'Save canceled'},handles.run_wipOutput)
    return
else
    data = guidata(handles.MainFigure);
    save([path '/' fileName],'data')
    setOutput({'Save done!'},handles.run_wipOutput)
end
