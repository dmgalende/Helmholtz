function plotBottom(handles)

data = guidata(handles.MainFigure);
if isfield(data,'plotBottom')
    delete(data.plotBottom)
end

if data.plotCondition
    hold on
    if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
        nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
        nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};
        data.plotBottom = plotSolution(data.mesh.(nameNodal),data.mesh.(nameCon),...
            data.bottom.value,data.mesh.referenceElement);
    elseif strcmp(data.computation,'NEFEM')
        data.plotBottom = plotSolution_NEFEM(data.mesh,data.bottom.value);
    end
    hold off
else
    data.plotBottom = [];
end

if ~get(handles.sel_bottomCheck,'Value')
    set(data.plotBottom,'Visible','off')
end

guidata(handles.MainFigure,data)
