function checkVisibility(nefemMeshFile)

home, close all

%----------------------- Options
plotInteriorElems = false;
npoints = 1000;
tolSelectionNode = 1e-1;
maxTangentAngle = 8*pi/180;

%----------------------- Load data and detect boundaries
addpath auxFunctions lib_nurbs

data = load(nefemMeshFile);
fields = fieldnames(data);
k = [];
kPos = size(data.T,1);
totalT = data.T;
for ifield = 1:numel(fields);
    iname = fields{ifield};
    if length(iname) > 3 && strcmp(iname(1:3),'Tb_')
        k = [k ifield];
        totalT = [totalT ; data.(iname)];
        kPos = [kPos kPos(end)+size(data.(iname),1)];
    end
end

%----------------------- Discretize nurbs in npoints
xcoordDis = linspace(-1,1,npoints); 
coordRefDis = [xcoordDis' -1*ones(size(xcoordDis'))]; %first face: [-1,1]x[-1,-1]

%----------------------- Check visibility for each boundary and each element
X = data.X;
nurbs = data.nurbs;
nElements = size(totalT,1);
elemProblemInfo = zeros(nElements,4);
iproblem = 0;
for boundaryIndex = k
    iname = fields{boundaryIndex};
    T = data.(iname);
    nOfElements = size(T,1);
    trimmedInfo = data.trimmedInfo.(iname);
    for elem = 1:nOfElements
        Te = T(elem,:);
        Xe = X(Te,:);
        
        [u,nurbsxy] = getVisibilityCondition();
        
        if u < npoints
            disp(['Problem in elem ' num2str(elem) ' and connectivity ' iname])
            iproblem = iproblem + 1;
            elemProblemInfo(iproblem,:) = [boundaryIndex elem nurbsxy];
        end
    end
end

%----------------------- Fix the mesh manually (if necesary)
if iproblem == 0
    disp(['The mesh has no visibility problems (checked with ' num2str(npoints) ' points)'])
    return
end
    
% GUI for fixing the nodes position by hand
disp('Initializing figure for fixing...')
MainFigure = figure(...
                    'Visible','off',...
                    'WindowButtonDownFcn', @MainFigureWindowButtonDownFcn,...
                    'WindowButtonUpFcn', @MainFigureWindowButtonUpFcn,...
                    'WindowButtonMotionFcn', @MainFigureWindowButtonMotionFcn);
axesHandle = axes();
axis equal

% Nodal connectivity
disp('Computing nodal connectivity...')
estimatedValence = 7;
nNodes = max(max(totalT));
totalN = zeros(nNodes,estimatedValence);
nn = ones(nNodes,1);
for ielem = 1:nElements
    Te = totalT(ielem,:);
    nn_Te = nn(Te);
    for kk = 1:3
        totalN(Te(kk),nn_Te(kk)) = ielem;
    end
    nn(Te) = nn(Te) + 1;
end
totalN(:,max(nn):end) = [];
clear nn

% Double ring of elements
disp('Computing double ring of boundary elements...')
doubleRing = zeros(estimatedValence*nElements,1);
aux1 = 1;
aux2 = 0;
for cont = k
    name = fields{cont};
    disp(['-- Boundary ' name(4:end)])
    extT = data.(name)(:,1:3);
    nOfElems = size(extT,1);
    for e = 1:nOfElems
        extNode = extT(e,3);
        intElem = totalN(extNode,logical(totalN(extNode,:)));
        aux2 = aux2 + length(intElem);
        doubleRing(aux1:aux2,1) = intElem;
        aux1 = aux2 + 1;
    end
end
doubleRing(aux2+1:end) = [];
doubleRing = unique(doubleRing);

% Plot interior elements
if plotInteriorElems
    patch('Faces',data.T(:,1:3),'Vertices',X,'FaceColor','none','EdgeAlpha',1);
else
    patch('Faces',totalT(doubleRing,1:3),'Vertices',X,'FaceColor','none','EdgeAlpha',1);
end

% Plot elements which have visibility problems
elemProblemInfo(iproblem+1:end,:) = [];
hElemProblem = zeros(iproblem,2);
hold on
for ielem = 1:iproblem
    boundaryIndex = elemProblemInfo(ielem,1);
    elem = elemProblemInfo(ielem,2);
    Te = data.(fields{boundaryIndex})(elem,1:3);
    hElemProblem(ielem,1) = patch('Faces',Te,'Vertices',X,...
                                  'FaceColor','r','EdgeAlpha',1,...
                                  'FaceAlpha',0.5);
    nurbPoint = elemProblemInfo(ielem,3:4);
    v3 = X(Te(3),:);
    hElemProblem(ielem,2) = plot([v3(1) nurbPoint(1)],[v3(2) nurbPoint(2)],...
                                 'g','lineWidth',2);
end

% Plot nurbs (boundary)
for inurb = 1:numel(nurbs)
    aNurb = nurbs(inurb);
    h = nurbsCurvePlot(aNurb,aNurb.iniParam,aNurb.endParam,3000);
    set(h,'color','k','lineWidth',2)
end
hold off

% Initialization control and auxiliar variables
userSelectionNode = false;
userSelectionElem = false;
selectedElem = 0;
fixedElems = false(iproblem,1);
fixedNodes = false(nNodes,1);
xy = [0 0];
hSelectedNode = [];
hSelectedLines = [];
rgbColorSelected = [1 1 0];
set(MainFigure,'Visible','on')
disp('Click on a red element in the figure, select the vertex and drag it to new position')
uiwait(MainFigure)

% Compute new nodes position
disp('Computing new nodes position...')
nodesVec = 1:nNodes;
nodesVec = nodesVec(fixedNodes);
elemVec = false(nElements,1);
referenceElement = createReferenceElement(data.elemInfo.type,data.elemInfo.nOfNodes,[]);
nodesCoord = referenceElement.NodesCoord;
for inode = nodesVec
    elems = totalN(inode,logical(totalN(inode,:)));
    for ielem = elems
        if elemVec(ielem), continue, end
        index = 1;
        while ielem > kPos(index), index = index + 1; end
        if index == 1 %interior elements
            Te = data.T(ielem,:);
            nodesAdapted = nodesCoord;
        else %boundary elements
            jelem = ielem - kPos(index-1);
            index = k(index-1);
            Te = data.(fields{index})(jelem,:);
            idNurbs = data.trimmedInfo.(fields{index})(jelem).idNurbs;
            trim = data.trimmedInfo.(fields{index})(jelem).trim;
            nodesAdapted = nefemInterp2DAdaptedNodesElement(nodesCoord,X(Te(1:3),:),...
                nurbs(idNurbs),trim(1),trim(2));
        end
        X(Te,:) = linearMapping(X(Te(1:3),:),nodesAdapted);
        elemVec(ielem) = true;
    end
end    

% Save new data
data.X = X;
save([nefemMeshFile(1:end-4) '_fixed'],'-struct','data')


%---------------------------------------------------------------
%----------------------- Callbacks and auxiliar nested functions
%---------------------------------------------------------------

    function MainFigureWindowButtonDownFcn(hObject,eventdata)
        if userSelectionElem
            cp = get(axesHandle,'currentPoint');
            xy = [cp(1,1) cp(1,2)];
            boundaryIndex = elemProblemInfo(selectedElem,1);
            elem = elemProblemInfo(selectedElem,2);
            Te = data.(fields{boundaryIndex})(elem,1:3);
            Xe = X(Te,:);
            if norm(Xe(3,:)-xy)/norm(Xe(3,:)) <= tolSelectionNode
                userSelectionNode = true;
                selectedNode = Te(3); 
                hold on, hSelectedNode = plot(xy(1),xy(2),'go',...
                    'markerFaceColor','g','markerSize',10);
                hold off
            else
                if fixedElems(selectedElem)
                    set(hElemProblem(selectedElem,1),'FaceColor','g')
                else
                    set(hElemProblem(selectedElem,1),'FaceColor','r')
                end
                userSelectionElem = false;
            end
        elseif strcmp(get(gco,'type'),'patch') 
            selectedElem = find(hElemProblem(:,1) == gco);
            if ~isempty(selectedElem)
                userSelectionElem = true;
                set(hElemProblem(selectedElem,1),'FaceColor',rgbColorSelected)
            end
        else
            disp('Click on a red element in the figure')
        end
    end

%------------

    function MainFigureWindowButtonUpFcn(hObject,eventdata)
        if userSelectionNode
            delete([hSelectedNode hSelectedLines])
            userSelectionElem = false;
            userSelectionNode = false;
            userAns = questdlg('Do you want to keep the new node position?',...
                'Sure?','Yes','No','No');
            if strcmp(userAns,'Yes')
                iname = fields{boundaryIndex};
                trimmedInfo = data.trimmedInfo.(iname);
                inode = data.(iname)(elem,3);
                X(inode,:) = xy;
                Xe = X(data.(iname)(elem,1:3),:);
                v3 = Xe(3,:);
                [u,nurbsxy] = getVisibilityCondition();
                delete(hElemProblem(selectedElem,:))
                hold on
                hElemProblem(selectedElem,1) = patch('Faces',Te,'Vertices',X,...
                    'EdgeAlpha',1,'FaceAlpha',0.5);
                if u < npoints
                    disp('The element has not been fixed')
                    hElemProblem(selectedElem,2) = plot([v3(1) nurbsxy(1)],...
                        [v3(2) nurbsxy(2)],'g','lineWidth',2);
                    set(hElemProblem(selectedElem,1),'FaceColor','r')
                    fixedElems(selectedElem) = false;
                    fixedNodes(inode) = false;
                else
                    disp('The element has been fixed')
                    hElemProblem(selectedElem,2) = plot(v3(1),v3(2),'g'); %avoid easily errors in line plot
                    set(hElemProblem(selectedElem,1),'FaceColor','g')
                    fixedElems(selectedElem) = true;
                    fixedNodes(inode) = true;
                end
                hold off
            else
                if fixedElems(selectedElem)
                    set(hElemProblem(selectedElem,1),'FaceColor','g')
                else
                    set(hElemProblem(selectedElem,1),'FaceColor','r')
                end
            end
            hSelectedNode = [];
            hSelectedLines = [];
        end
    end

%------------

    function MainFigureWindowButtonMotionFcn(hObject,eventdata)
        if userSelectionNode
            cp = get(axesHandle,'currentPoint');
            xy = [cp(1,1) cp(1,2)];
            lpoints = [Xe(1,:) ; xy ; Xe(2,:)];
            delete([hSelectedNode hSelectedLines])
            hold on, hSelectedNode = plot(xy(1),xy(2),'go',...
                'markerFaceColor','g','markerSize',10);
            hSelectedLines = plot(lpoints(:,1),lpoints(:,2),'y',...
                'lineWidth',2);
            hold off
        end
    end

%------------

    function [u,nurbsFinalPoint] = getVisibilityCondition()
        
        % NEFEM information
        u1 = trimmedInfo(elem).trim(1);
        u2 = trimmedInfo(elem).trim(2);
        idNurbs = trimmedInfo(elem).idNurbs;
        aNurbs = nurbs(idNurbs);
        vertCoord = Xe(1:3,:);
        v3 = vertCoord(3,:);
        
        % Adapted reference element nodes (first face) to nurbs size
        nodesAdapted = nefemInterp2DAdaptedNodesElement...
            (coordRefDis,vertCoord,aNurbs,u1,u2);
        
        % Map nodes to real space
        nurbsxy = linearMapping(vertCoord,nodesAdapted);
        
        % Check angle between lines from vertex v3
        prev_angle = 0;
        l1 = v3 - nurbsxy(1,:);
        l1 = l1/norm(l1);
        l2 = v3 - nurbsxy(2,:);
        new_angle = acos(l1*l2'/norm(l2));
        u = 2;
        while new_angle > prev_angle && u < npoints
            prev_angle = new_angle;
            l2 = v3 - nurbsxy(u+1,:);
            new_angle = acos(l1*l2'/norm(l2));
            u = u + 1;
        end
        
        nurbsFinalPoint = nurbsxy(u,:);
        
        if u == npoints
            v2 = vertCoord(2,:);
            v1 = vertCoord(1,:);
            l1 = v3 - v1; l2 = v3 - v2;
            l1 = l1/norm(l1); l2 = l2/norm(l2);
            auxL1 = nurbsCurveDerivPoint(aNurbs,u1);
            auxL2 = nurbsCurveDerivPoint(aNurbs,u2);
            L1 = auxL1(1:2)/norm(auxL1); L2 = auxL2(1:2)/norm(auxL2);
            prod1 = abs(l1*L1'); prod2 = abs(l2*L2');
            if any([prod1,prod2] > cos(maxTangentAngle))
                u = u - 1;
            end
        end
    end

%------------

end

        
        
        
        
    
