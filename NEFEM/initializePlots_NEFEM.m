function initializePlots_NEFEM(handles)

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
nurbs = data.mesh.nurbs;
nOfBoundaries = numel(data.mesh.boundaryNames);
data.plotBoundary = cell(1,nOfBoundaries);
for i = 1:nOfBoundaries
    item = data.mesh.boundaryIndex(i);
    iboundary = data.mesh.fieldNames{item};
    iboundaryTrimInfo = data.mesh.trimmedInfo.(iboundary);
    nOfElements = size(data.mesh.(iboundary),1);
    data.plotBoundary{i} = zeros(nOfElements,1);
    hold on
    for j = 1:nOfElements
        u1 = iboundaryTrimInfo(j).trim(1);
        u2 = iboundaryTrimInfo(j).trim(2);
        idNurbs = iboundaryTrimInfo(j).idNurbs;
        aNurbs = nurbs(idNurbs);
        h = nurbsCurvePlot(aNurbs,u1,u2,50);
        set(h,'lineWidth',3,'color','r')
        data.plotBoundary{i}(j) = h;
    end
    hold off
end
axis equal

% Mesh and submeshes
if data.plotCondition
    hold on
    
    % Mesh
    data.plotMesh = plotMesh_NEFEM(data.mesh,30);

    % Submeshes (treated as standard FEM plots)
    nOfSubmeshes = numel(data.mesh.submeshNames);
    for i = 1:nOfSubmeshes
        item = data.mesh.submeshIndex(i);
        imesh = data.mesh.fieldNames{item};
        data.plotSubmesh(i) = plotMesh(data.mesh.X,data.mesh.(imesh),...
            data.mesh.elemInfo.faceNodes);
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
