function elemType = getReferenceElement(handles)

data = guidata(handles.MainFigure);

if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
    nameElem = data.mesh.fieldNames{data.mesh.indexElemPosCon(1)};
    nameType = data.mesh.elemFieldNames{data.mesh.indexTypNodFac1Fac2(1)};
    nameNodes = data.mesh.elemFieldNames{data.mesh.indexTypNodFac1Fac2(2)};
elseif strcmp(data.computation,'NEFEM')
    nameElem = 'elemInfo';
    nameType = 'type';
    nameNodes = 'nOfNodes';
end
if data.mesh.(nameElem).(nameType)
    elemString = 'Tri ';
else
    elemString = 'Qua ';
end
elemType = [elemString 'with ' num2str(data.mesh.(nameElem).(nameNodes))...
    ' nodes'];
data.mesh.referenceElement = createReferenceElement(data.mesh.(nameElem).(nameType),...
    data.mesh.(nameElem).(nameNodes),handles.run_wipOutput);

guidata(handles.MainFigure,data);