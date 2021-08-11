function initializePlots(handles)

data = guidata(handles.MainFigure);
cla(handles.axesHandle,'reset')
if isfield(data,'plotBoundary')
    data = rmfield(data,{'plotBoundary' 'plotMesh'});
    if isfield(data,'plotSubmesh')
        data = rmfield(data,'plotSubmesh');
    end
    if isfield(data,'plotBottom')
        data = rmfield(data,'plotBottom');
    end
end

% Boundary
nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
nOfBoundaries = numel(data.mesh.boundaryNames);
data.plotBoundary = cell(1,nOfBoundaries);
for i=1:nOfBoundaries
    item = data.mesh.boundaryIndex(i);
    iboundary = data.mesh.fieldNames{item};
    nOfElements = size(data.mesh.(iboundary),1);
    data.plotBoundary{i} = zeros(nOfElements,1);
    hold on
    for j=1:nOfElements
        Tf = data.mesh.(iboundary)(j,:);
        Xf = data.mesh.(nameNodal)(Tf,:);
        data.plotBoundary{i}(j) = plot(Xf(:,1),Xf(:,2),'r','LineWidth',3);
    end
    hold off
end
axis equal

% Mesh and submeshes
nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};
nameElem = data.mesh.fieldNames{data.mesh.indexElemPosCon(1)};
nameFac2d = data.mesh.elemFieldNames{data.mesh.indexTypNodFac1Fac2(4)};
nOfSubmeshes = numel(data.mesh.submeshNames);
if data.plotCondition
    hold on

    % Mesh
    data.plotMesh = plotMesh(data.mesh.(nameNodal),data.mesh.(nameCon),...
        data.mesh.(nameElem).(nameFac2d));

    % Submeshes
    data.plotSubmesh = zeros(nOfSubmeshes,1);
    for i = 1:nOfSubmeshes
        item = data.mesh.submeshIndex(i);
        imesh = data.mesh.fieldNames{item};
        data.plotSubmesh(i) = plotMesh(data.mesh.(nameNodal),data.mesh.(imesh),...
            data.mesh.(nameElem).(nameFac2d));
        set(data.plotSubmesh(i),'FaceColor',[0 0 1],'Visible','off')
    end
    hold off
else
    data.plotMesh = [];
    data.plotSubmesh = [];
end

% Visibility 
if ~get(handles.sel_meshCheck,'Value')
    set(data.plotMesh,'Visible','off')
end

if ~get(handles.sel_axesCheck,'Value')
    set(gca,'Visible','off')
end

guidata(handles.MainFigure,data);