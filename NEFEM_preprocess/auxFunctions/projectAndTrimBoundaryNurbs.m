function [mesh,trimmedInfo] = projectAndTrimBoundaryNurbs(mesh,meshFields,k)

global nurbs projectNodes referenceElement adaptedNodesToNurbsBoundary

TOLproj = 5e-5;
X = mesh.X;
nOfTotalNodes = size(X,1);
intT = mesh.T;
intTLinear = intT(:,1:3);
coordRef = referenceElement.NodesCoord;
faceVertices = computeFaceVertices2D();
refnOfNodes = mesh.elemInfo.nOfNodes;
mapFace2 = load(['rotTriangle_' num2str(refnOfNodes) 'nodes_face2']);
mapFace3 = load(['rotTriangle_' num2str(refnOfNodes) 'nodes_face3']);
for i = k
    disp('-----------------')
    disp(['Connectivity ' meshFields{i}])
    disp('-----------------')
    T = mesh.(meshFields{i});
    
    %Boundary nodes
    nodesCurvedBoundary = computeNurbsBoundaryNodes(X,T,faceVertices);
    
    %Project the boundary nodes
    if projectNodes
        disp('Projecting nodes...')
        nOfNodes = size(nodesCurvedBoundary,2);
        for kNode = 1:nOfNodes
            iNode = nodesCurvedBoundary(1,kNode);
            iNurbs = nodesCurvedBoundary(2,kNode);
            xNode = X(iNode,:);
            aNurbs = nurbs(iNurbs);
            pt1 = nurbsCurvePoint(aNurbs, aNurbs.iniParam);
            pt2 = nurbsCurvePoint(aNurbs, aNurbs.endParam);
            if norm(xNode - pt1(1:2))<TOLproj
                ptProj = pt1(1:2);
            elseif norm(xNode - pt2(1:2))<TOLproj
                ptProj = pt2(1:2);
            else
                uIni = initialGuessProjection(aNurbs, xNode);
                paramNode = fzero(@(u) fForZero(u, xNode, aNurbs), uIni);
                ptProj = nurbsCurvePoint(aNurbs, paramNode);
            end
            X(iNode,:) = ptProj(1:2);
        end
    end
    
    %Trim
    nOfElements = size(T,1);
    disp(['TRIM PROCESS: ' num2str(nOfElements) ' elements'])
    nodeCheck = false(nOfTotalNodes,1);
    aux = cell(nOfElements,1);
    aux2 = aux;
    aux(:) = {0};
    aux2(:) = {[0 0]};
    trimStruct = struct('idNurbs',aux(:),'trim',aux2(:)); %initialize structure
    aux = nodesCurvedBoundary(2,:);
    for elem = 1:nOfElements
        fprintf('-- -- -- element %i\n',elem);
        
        face = T(elem,end);
        faceNodesPos = faceVertices(face,:);
        nodes = T(elem,faceNodesPos);
        Xf = X(nodes,1);
        Yf = X(nodes,2);
        nurb1 = aux(nodesCurvedBoundary(1,:) == nodes(1));
        nurb2 = aux(nodesCurvedBoundary(1,:) == nodes(2));
        if ~isscalar(nurb1) && ~isscalar(nurb2)
            condNurb = (nurb1(1) == nurb2) | (nurb1(2) == nurb2);
            boundaryNurbs = nurb2(condNurb);
        elseif ~isscalar(nurb1)
            condNurb = nurb1 == nurb2;
            boundaryNurbs = nurb1(condNurb);
        else
            condNurb = nurb1 == nurb2;
            boundaryNurbs = nurb2(condNurb);
        end
        aNurbs = nurbs(boundaryNurbs);
        u1 = nurbsCurvePointInversion([Xf(1),Yf(1)],aNurbs);
        if u1 == -1 && ~projectNodes
            uIni = initialGuessProjection(aNurbs, [Xf(1),Yf(1)]);
            %u1 = fzero(@(u) fForZero(u, [Xf(1),Yf(1)], aNurbs), uIni);
            u1 = uIni;
        end   
        u2 = nurbsCurvePointInversion([Xf(2),Yf(2)],aNurbs);
        if u2 == -1 && ~projectNodes
            uIni = initialGuessProjection(aNurbs, [Xf(2),Yf(2)]);
            %u2 = fzero(@(u) fForZero(u, [Xf(2),Yf(2)], aNurbs), uIni);
            u2 = uIni;
        end
        if u1 == u2
            error(['nurbsCurvePointInversion takes the same parameter, consider'...
                   'increasing the tolerance value'])
        end
        trimStruct(elem).idNurbs = boundaryNurbs;
        trimStruct(elem).trim(1) = u1;
        trimStruct(elem).trim(2) = u2;
        
        %Renumber nodes: first face always curved
        if face == 2
            T(elem,1:end-1) = T(elem,mapFace2.rotationMap);
        elseif face == 3
            T(elem,1:end-1) = T(elem,mapFace3.rotationMap);
        end

        if adaptedNodesToNurbsBoundary
            %Map the reference element nodes to real space using the linear mapping
            p1 = nurbsCurvePoint(aNurbs,u1);
            p2 = nurbsCurvePoint(aNurbs,u2);
            vertCoord = [p1(1:2) ; p2(1:2) ; X(T(elem,3),:)];

            %Adapt nodes to nurb side before mapping
            nodesAdapted = nefemInterp2DAdaptedNodesElement(coordRef,vertCoord,aNurbs,u1,u2);

            %Map to real space and store in nodal matrix
            %         X(T(elem,:),:) = linearMapping(vertCoord,coordRef);
            X(T(elem,1:end-1),:) = linearMapping(vertCoord,nodesAdapted);

            %Recompute nodes position of those elements conected with elem
            p12 = [p1 ; p2];
            for cont = 1:2
                node = nodes(cont);
                if ~nodeCheck(node)
                    [intElem,intPos] = find(intTLinear == node);
                    ind = 1;
                    for ielem = intElem' %transpose is important!
                        ipos = intPos(ind);
                        vertexElem = X(intTLinear(ielem,:),:);
                        vertexElem(ipos,:) = p12(cont,1:2);
                        X(intT(ielem,:),:) = linearMapping(vertexElem,coordRef);
                        ind = ind + 1;
                    end
                    nodeCheck(node) = true;
                end
            end
        end
    end
    
    trimmedInfo.(meshFields{i}) = trimStruct;
    T(:,end) = []; %The face info is no more necessary
    mesh.(meshFields{i}) = T;
end
mesh.X = X;

    
    
    