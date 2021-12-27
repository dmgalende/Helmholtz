function fixBottom_Callback(hObject,eventdata,handles)

switch get(hObject,'state')
    case 'on'
        data = guidata(handles.MainFigure);

        %Check PML condition
        PML = [];
        nOfBoundaries = numel(data.mesh.boundaryNames);
        totalPML = false(1,nOfBoundaries);
        boundaries = 1:nOfBoundaries;
        for iboundary = 1:nOfBoundaries
            if strcmp(data.PML{5,iboundary},'on')
                PML = [PML iboundary];
            end
        end
        if isempty(PML)
            setOutput({'??? No PML boundary has been assigned yet'},handles.run_wipOutput)
            set(hObject,'state','off')
            return
        end

        %Turn off visibility of mesh, PML nodes and bottom plots
        handles2BnoVisible = [data.plotMesh data.plotBottom...
            data.plotPMLnodes{1,:} data.plotPMLnodes{2,:}];
        if ~isempty(data.mesh.submeshIndex) %Add submesh if it exists
            handles2BnoVisible = [handles2BnoVisible data.plotSubmesh'];
        end
        handles2Bvisible = [];
        for iboundary = 1:numel(data.mesh.boundaryNames)
            auxHandle = data.plotPMLarea{1,iboundary};
            handles2Bvisible = [handles2Bvisible auxHandle'];
        end
        prevVisHandles = get([handles2BnoVisible handles2Bvisible],'visible');
        newColor = [0 1 0]; %Color of PML selection (green)
        set(handles2BnoVisible,'visible','off')
        set(handles2Bvisible,'visible','on')

        %Waiting to user's action click
        setOutput({'Click in PML on the figure'},handles.run_wipOutput)
        set(handles.MainFigure,'WindowButtonDownFcn',@axesClicking)
        uiwait(handles.MainFigure)

        %Read selected PMLs and try to do the smoothing process
        set(handles.MainFigure,'WindowButtonDownFcn',[])
        totalPML = boundaries(totalPML);
        if isempty(totalPML)
            setOutput({'No PMLs selected. Process canceled'},...
                handles.run_wipOutput)
            return2previousPlots %Return to pre-process plot state
        else
            msje = {'hconst = constant bottom on PML area'
                'alpha = maximum distance of smoothing process'
                'Enter the 2x1 matrix [hconst ; alpha] for each selected PML'};
            setOutput(msje,handles.run_wipOutput)
            nOfPMLs = length(totalPML);
            inputfun %Ask user and create the parameters array (cell matrix)
            if isempty(parameters)
                setOutput({'Process canceled by user'},handles.run_wipOutput)
                return2previousPlots %Return to pre-process plot state
                return
            end
            computeSmoothBottom %run smooth process and create the bottom array newBottom
            if data.plotCondition
                plotNewBottom(nameNodal,T_delaunay,newBottom)
                answer = questdlg('Do you want to keep the new bottom data?',...
                    'Process finished','Yes','No','Yes');
                if strcmp(answer,'Yes')
                    data.bottom.value = newBottom;
                    delete(data.plotBottom)
                    data.plotBottom = plotNewBottomHandle;
                    aux = numel(data.plotMesh);
                    handles2BnoVisible(aux+1:2*aux+1) = plotNewBottomHandle;  %New bottom patch and colorbar
                    if isfield(data.bottom,'ccg')
                        data.bottom = rmfield(data.bottom,{'waveNumber' 'ccg'});
                    end
                    setOutput({'New bottom data has been stored'},handles.run_wipOutput)
                else
                    delete(plotNewBottomHandle)
                    setOutput({'No data has been stored'},handles.run_wipOutput)
                end
            else
                data.bottom.value = newBottom;
                if isfield(data.bottom,'ccg')
                    data.bottom = rmfield(data.bottom,{'waveNumber' 'ccg'});
                end
                setOutput({'New bottom data has been stored'},handles.run_wipOutput)
            end

            guidata(handles.MainFigure,data);

            return2previousPlots %Return to pre-process plot state
        end
    case 'off'
        setOutput({'PML selection disabled. Click on figure to continue'},...
            handles.run_wipOutput)
end


%_________________
    function axesClicking(hObject,eventdata)

        cp = get(handles.axesHandle,'currentPoint');
        xy = cp(1,1:2);
        selectedPML = findPML(xy);
        moreClicks = strcmp(get(handles.fixBottom,'state'),'on');
        if ~moreClicks
            uiresume(handles.MainFigure)
        elseif ~selectedPML && moreClicks
            setOutput({['PML not selected. Try again or disable the tool by '...
                'cliking on the icon']},handles.run_wipOutput)
        elseif totalPML(selectedPML)
            setOutput({['PML not selected on boundary "' data.mesh.boundaryNames{selectedPML}...
                '". Click on more PML or disable the tool by cliking on the icon']},...
                handles.run_wipOutput)
            set(data.plotPMLarea{1,selectedPML},'faceColor','r')
            totalPML(selectedPML) = false;
        else
            setOutput({['PML selected on boundary "' data.mesh.boundaryNames{selectedPML}...
                '". Click on more PML or disable the tool by cliking on the icon']},...
                handles.run_wipOutput)
            set(data.plotPMLarea{1,selectedPML},'faceColor',newColor)
            totalPML(selectedPML) = true;
        end

    end
%_________________


%_________________
    function selectedPML = findPML(xy)

        selectedPML = 0;
        for iPML = PML
            sigma = coefPML(xy,...
                1,... %do xy(1)                
                data.PML{2,iPML},...
                1,... %value ignored
                1,... %value ignored
                data.BC.parameters{iPML}{4}(3));
            if any(sigma ~= 0)
                selectedPML = iPML;
                break
            end
        end

    end
%_________________


%_________________
    function inputfun

        PMLnames = data.mesh.boundaryNames(totalPML);
        parameters = inputdlg(PMLnames,'Parameters',2);

        %Checking cancel
        if isempty(parameters)
            return
        end

        %Checking parameters
        for j = 1:nOfPMLs
            jparam = parameters{j};
            if ~isempty(jparam)
                jparam = str2num(jparam);
                if any(size(jparam) ~= [2 1])
                    setOutput({'Invalid input. Remind that it has to be a 2x1 numerical matrix'},...
                        handles.run_wipOutput)
                    inputfun
                    break
                else
                    parameters{j} = jparam;
                end
            else
                setOutput({'Invalid input. Remind that it has to be a 2x1 numerical matrix'},...
                    handles.run_wipOutput)
                inputfun
                break
            end
        end

    end
%_________________


%_________________
    function return2previousPlots

        sizeNoVisible = length(handles2BnoVisible);
        sizeVisible = length(handles2Bvisible);

        %Mesh, bottom, submehes and PML nodes plots
        for j = 1:sizeNoVisible
            set(handles2BnoVisible(j),'visible',prevVisHandles{j})
        end

        %PML area
        for j = 1:sizeVisible
            set(handles2Bvisible(j),'visible',prevVisHandles{j+sizeNoVisible},...
                'faceColor','r')
        end

    end
%_________________


%_________________
    function computeSmoothBottom

        %Mesh info
        if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
            nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
            nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};
            conec = data.mesh.(nameCon);
        elseif strcmp(data.computation,'NEFEM')
            nameNodal = 'X';
            nOfBoundaries = numel(data.mesh.boundaryNames);
            conec = data.mesh.T;
            for iboundaryNEFEM = 1:nOfBoundaries
                inameCon = data.mesh.fieldNames{data.mesh.boundaryIndex(iboundaryNEFEM)};
                conec = [conec ; data.mesh.(inameCon)];
            end
        end
        X = data.mesh.(nameNodal);
        T_delaunay = plotSolution(X,conec,0,data.mesh.referenceElement,true); %Delaunay mesh (linear)
        N = getNodalConnectivity(T_delaunay); %Nodal connectivity

        %PML info
        PMLinfo = cell(1,4,nOfPMLs); %Coordinates of each cartesian PML boundary
        PMLlon = zeros(nOfPMLs,1); %Length of each PML area
        for i = 1:nOfPMLs
            iPML = data.PML{2,totalPML(i)};
            PMLinfo(:,:,i) = iPML;
            lon = data.BC.parameters{totalPML(i)}{4}(3);
            PMLlon(i) = lon;
        end

        %---------------------
        %INCOMING FRONT SEARCH
        %---------------------

        setOutput({'Beginnig incoming front search...'},handles.run_wipOutput)

        %PML nodes
        sizeX = size(X,1);
        numbering = 1:sizeX; %vector of mesh numbering
        [nodesPML,T,newBottom] = getPMLnodesAndInitialBottom(T_delaunay,N,sizeX,numbering); %nodesPML is a boolean array
        %T is the delaunay linear mesh with PML nodes = -1
        
        %Initialization
        front = nodesPML;
        marked = true;
        nodes2modify = false(sizeX,1);
        distances2store = zeros(sizeX,2);
        plotFrontHandle = [];
        k = 1;

        %Algorithm (boolean)
        while any(marked)

            %Compute new front and marked elements
            [neighbors,T] = findFrontNeighbors(front,sizeX,T,N,numbering);
            [marked,distances] = findMinDistances(neighbors,PMLinfo,PMLlon,X,numbering);
            plotFrontHandle = plotFront(X,neighbors,plotFrontHandle);
            setOutput({['Front ' num2str(k) ': ' num2str(length(X(marked,1)))...
                ' new nodes are considered']},handles.run_wipOutput);

            %Update front and data
            front = neighbors & marked;
            nodes2modify = nodes2modify | marked;
            distances2store = distances2store + distances;
            k = k + 1;
        end

        %New bottom
        setOutput({'Computing new bottom...'},handles.run_wipOutput)
        nodesNewBottom = numbering(nodes2modify);
        for inode = nodesNewBottom
            dis = distances2store(inode,1);
            hconst = parameters{distances2store(inode,2)}(1);
            alpha = parameters{distances2store(inode,2)}(2);

            newBottom(inode) = newBottom(inode)*dis/alpha +...
                hconst*(alpha - dis)/alpha;
        end
        delete(plotFrontHandle)

    end
%_________________


%_________________
    function [nodesPMLi,T_mod,newBottom] = getPMLnodesAndInitialBottom(T,N,sizeX,numbering)

        %Initialize newBottom
        newBottom = data.bottom.value;

        %PML nodes and constant bottom
        nodesPML = false(sizeX,1);
        for i = 1:nOfPMLs
            nodesPMLaux = false(sizeX,1);
            param = data.BC.parameters{totalPML(i)}{4}(4);
            hconst = parameters{i}(1);
            if param > 1
                nameSubCon = data.mesh.fieldNames{data.mesh.submeshIndex(param-1)};
                nodesPMLaux(getNodes(data.mesh.(nameSubCon))) = true;
            else
                nodesPMLaux(data.PML{4,totalPML(i)}{1}) = true;
                nodesPMLaux(data.PML{4,totalPML(i)}{2}) = true;
            end
            newBottom(nodesPMLaux) = hconst;
            nodesPML = nodesPML | nodesPMLaux;
        end

        %Conectivity without PML elements
        nodesPMLnum = numbering(nodesPML);
        nOfElements = size(T,1);
        nOfTimesMarked = zeros(nOfElements,4);
        for inode = nodesPMLnum
            elem = N(inode,logical(N(inode,:)));
            for j = 1:length(elem)
                jelem = elem(j);
                jpos = find(T(jelem,:) == inode);
                nOfTimesMarked(jelem,1) = nOfTimesMarked(jelem,1) + 1;
                nOfTimesMarked(jelem,jpos+1) = 1;
            end
        end
        beRemove = nOfTimesMarked(:,1) == 3; %Only elements with all PML nodes
        T_mod = T;
        T_mod(beRemove,:) = -1;

        %PML interfase's nodes
        beConsidered = ~beRemove & (nOfTimesMarked(:,1) > 0);
        elemNum = 1:nOfElements;
        elemConsidered = elemNum(beConsidered);
        nodesPMLi = false(sizeX,1);
        for ielem = elemConsidered
            ielemCond = nOfTimesMarked(ielem,2:end) == 1;
            nodesPMLi(T(ielem,ielemCond)) = true;
        end

    end
%_________________


%_________________
    function [neighbors,T] = findFrontNeighbors(front,sizeX,T,N,numbering)

        %Initialization
        neighbors = false(sizeX,1);

        %Find neighbors
        frontNodes = numbering(front);
        for inode = frontNodes
            elem = N(inode,logical(N(inode,:)));
            for j = 1:length(elem)
                jelem = elem(j);
                ipos = (T(jelem,:) > 0) & (T(jelem,:) ~= inode);
                neighbors(T(jelem,ipos)) = true;
            end
            T(elem,:) = -1;
        end
        neighbors = xor(neighbors,neighbors & front); %Exclusive neighbors

    end
%_________________


%_________________
    function [marked,finalDistances] = findMinDistances(neighbors,PMLinfo,...
            PMLlon,X,numbering)

        %Initialization
        sizeX = size(X,1);
        marked = false(sizeX,1);
        distances = NaN*ones(sizeX,1); %initialize distances at infinity
        finalDistances = zeros(sizeX,2);
        coord2check = [1 2 1 2];
        coef = [1 -1 -1 1];
        coef2 = [1 -1 1 -1];
        neighNodes = numbering(neighbors);

        %For each neighbor
        for inode = neighNodes
            xy = X(inode,:);

            %For each PML
            for iPML = 1:nOfPMLs
                alpha = parameters{iPML}(2);
                flag = false;

                %For each type of cartesian boundary
                for j = 1:4
                    jtype = PMLinfo{:,j,iPML};
                    jcoef = coef(j);
                    jcoef2 = coef2(j);
                    jsize = size(jtype,1);
                    jcoord = coord2check(j);
                    jlon = jcoef*PMLlon(iPML);

                    %For each straight line of the cartesian boundary
                    for k = 1:jsize

                        %Coordinates of the current line
                        ktype = jtype(k,:);
                        endPML = ktype(3) + jlon;

                        %Define the eight candidate's areas
                        area = defineAreasFromPMLline(ktype,endPML,alpha,j);
                        x = xy(jcoord);
                        y = xy(3 - jcoord);
                        xyaux = [x y];

                        %Compute minimum distance from rectangular boxes
                        for box = 1:2:8
                            x1 = area(box,1); x2 = area(box,2);
                            y1 = area(box,3); y2 = area(box,4);
                            if x >= x1 && x <= x2 && jcoef2*y >= jcoef2*y1 &&...
                                    jcoef2*y <= jcoef2*y2
                                point = area(box,5:6); %distance reference point
                                coord = area(box,7);
                                dis = abs(xyaux(coord) - point(coord));
                                prevDis = distances(inode);
                                distances(inode) = min(dis,prevDis);
                                flag = true;
                                break
                            end
                        end

                        %Compute minimum distance from cercle boxes
                        for box = 2:2:8
                            point = area(box,5:6); %center point
                            if (x-point(1))^2 + (y-point(2))^2 - alpha^2 <= 0
                                dis = norm(xyaux - point);
                                prevDis = distances(inode);
                                distances(inode) = min(dis,prevDis);
                                flag = true;
                            end
                        end

                        %Save distance
                        if flag
                            finalDistances(inode,:) = [distances(inode) iPML];
                            marked(inode) = true;
                        end
                    end
                end
                if flag, break, end
            end
        end

    end
%_________________


%_________________
    function area = defineAreasFromPMLline(x,yf,a,j)

        %Initialize
        area = zeros(8,7);

        %Change in coordinates to define a unique reference box
        if any(j == [3 4])
            aux = x(3);
            x(3) = yf;
            yf = aux;
        end

        %Change in length to define a unique reference box
        if any(j == [2 4])
            a = -a;
        end

        %Define reference box
        area(1,:) = [x(1) x(2) x(3)-a x(3) x(2) x(3) 2];
        area(2,:) = [x(2) x(2)+a x(3)-a x(3) x(2) x(3) 0];
        area(3,:) = [x(2) x(2)+a x(3) yf x(2) yf 1];
        area(4,:) = [x(2) x(2)+a yf yf+a x(2) yf 0];
        area(5,:) = [x(1) x(2) yf yf+a x(1) yf 2];
        area(6,:) = [x(1)-a x(1) yf yf+a x(1) yf 0];
        area(7,:) = [x(1)-a x(1) x(3) yf x(1) x(3) 1];
        area(8,:) = [x(1)-a x(1) x(3)-a x(3) x(1) x(3) 0];

    end
%_________________


%_________________
    function N = getNodalConnectivity(T_delaunay)

        nNodes = max(max(T_delaunay));
        nNodesElem = size(T_delaunay,2);
        N = zeros(nNodes,10);
        nn = ones(nNodes,1);
        for ielem = 1:size(T_delaunay,1)
            Te = T_delaunay(ielem,:);
            nn_Te = nn(Te);
            for kk = 1:nNodesElem
                N(Te(kk),nn_Te(kk)) = ielem;
            end
            nn(Te) = nn(Te) + 1;
        end
        N(:,max(nn):end) = [];

    end
%_________________


%_________________
    function plotFrontHandle = plotFront(X,neighbors,plotFrontHandle)

        delete(plotFrontHandle)
        hold on
        plotFrontHandle = plot(X(neighbors,1),X(neighbors,2),...
            'ro','markerSize',3,'markerFaceColor','r');
        hold off
        set(plotFrontHandle,'Parent',handles.axesHandle);

    end
%_________________


%_________________
    function plotNewBottom(nameNodal,T_delaunay,newBottom)

        hold on
        if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
            plotNewBottomHandle = plotSolution(data.mesh.(nameNodal),T_delaunay,newBottom);
        elseif strcmp(data.computation,'NEFEM')
            plotNewBottomHandle = plotSolution_NEFEM(data.mesh,newBottom,T_delaunay,-1);
        end

        hold off

        %Visibility
        set(handles2Bvisible,'visible','off')
        set(plotNewBottomHandle,'visible','on')

    end
%_________________


end






