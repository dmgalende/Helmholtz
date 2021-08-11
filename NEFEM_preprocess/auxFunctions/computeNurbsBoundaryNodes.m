function nurbsNodes = computeNurbsBoundaryNodes(X,T,faceVertices)

% nurbsnodes: only for the element vertices. The element has to belong to
% one boundary at the same time as maximum.

global nurbs referenceElementTri3

nExtElem = size(T,1);
exteriorNodes = zeros(nExtElem*2,1);
ini = 1;
fin = 2;
for iExtElem = 1:nExtElem
    iFace = T(iExtElem,end);
    faceNodes = faceVertices(iFace,:);
    exteriorNodes(ini:fin) = T(iExtElem,faceNodes);
    ini = fin + 1;
    fin = fin + 2;
end
exteriorNodes = unique(exteriorNodes);
ElemNodesInfo = createElemNodesInfo(T,exteriorNodes,faceVertices);

aux = [1 2];
x = X(exteriorNodes,1);
y = X(exteriorNodes,2);
tol = 1e-5;
nurbsNodes = [];
for iNurbs = 1:length(nurbs)
    
    % Check nurb has the initial and ending point on the boundary
    aNurb = nurbs(iNurbs);
    iniPoint = nurbsCurvePoint(aNurb,aNurb.iniParam);
    endPoint = nurbsCurvePoint(aNurb,aNurb.endParam);
    inicond1 = iniPoint(1) < x + tol & iniPoint(1) > x - tol;
    inicond2 = iniPoint(2) < y + tol & iniPoint(2) > y - tol;
    endcond1 = endPoint(1) < x + tol & endPoint(1) > x - tol;
    endcond2 = endPoint(2) < y + tol & endPoint(2) > y - tol;
    iniAndEndNode = exteriorNodes((inicond1 & inicond2) | (endcond1 & endcond2));
    nOfDetectedNodes = length(iniAndEndNode);
    if nOfDetectedNodes > 2
        error('POSSIBLE DOUBLED NODES DETECTED! CHECK THE MESH')
    elseif nOfDetectedNodes < 2
        continue
    end
    
    % Select direction for ordering nodes
    node = iniAndEndNode(1); %random node: initial or ending
    elemPos = 1; %random element
    u = 1;
    auxExteriorNodes = zeros(size(exteriorNodes));
    auxExteriorNodes(u) = node;
    elem = ElemNodesInfo(node,elemPos);
    iFace = T(elem,end);
    faceNodes = faceVertices(iFace,:);
    nodes = T(elem,faceNodes);
    inode = nodes(nodes ~= node);
    paramCheck = nurbsCurvePointInversion(X(inode,:),aNurb);
    if paramCheck == -1 %Make sure that direction for ordering nodes is the correct one
        elemPos = 2;
        elem = ElemNodesInfo(node,elemPos);
        if elem == 0 %Fix a possible bug when selecting discontinuous EZ4U boundaries
            continue
        end
        iFace = T(elem,end);
        faceNodes = faceVertices(iFace,:);
        nodes = T(elem,faceNodes);
        inode = nodes(nodes ~= node);
    end
    
    disp(['Nurb ' num2str(iNurbs) '...'])
    
    % Ordering nodes
    endNode = iniAndEndNode(2);
    while inode ~= endNode
        u = u + 1;
        auxExteriorNodes(u) = inode;
        elemPos = aux(ElemNodesInfo(inode,1:2) ~= elem);
        elem = ElemNodesInfo(inode,elemPos);
        iFace = T(elem,end);
        faceNodes = faceVertices(iFace,:);
        nodes = T(elem,faceNodes);
        inode = nodes(nodes ~= inode);
    end
    
    % Get the nurbs nodes info
    auxnodes = [auxExteriorNodes(1:u) ; endNode];
    nurbsNodes = [nurbsNodes  [auxnodes' ; iNurbs*ones(size(auxnodes'))]];
    
    % Make sure that nurb is defined with the same direction as the mesh (same tangent vector)
    tmesh = referenceElementTri3.N1dxi(1,:)*X(nodes,:); %using last element of ordering nodes process
    tparam = nurbsCurvePointInversion(X(endNode,:),aNurb);
    tnurb = nurbsCurveDerivPoint(aNurb,tparam);
    tmesh = tmesh/norm(tmesh);
    tnurb = tnurb(1:2)/norm(tnurb(1:2));
    cond = acos(tmesh*tnurb');
    if cond == pi/2 %Undefined case
        iniNode = nodes(nodes ~= endNode);
        tparam2 = initialGuessProjection(aNurb,X(iniNode,:));
        h = (tparam2 - tparam)/20;
        cont = 1;
        while cond == pi/2 && cont <= 21
            disp(['Angle between elem ' num2str(elem) ' and nurb ' num2str(iNurbs)...
                ' is pi/2. Taking another point to check'])
            tnurb = nurbsCurveDerivPoint(aNurb,tparam+cont*h);
            tnurb = tnurb(1:2)/norm(tnurb(1:2));
            cond = acos(tmesh*tnurb');
            cont = cont + 1;
        end
    end
    if cond > pi/2 %Change nurb definition sense
        disp('Changing definition sense...')
        order = size(aNurb.Pw,1);
        auxPw = zeros(size(aNurb.Pw));
        for ipoint = 1:order
            auxPw(order-ipoint+1,:) = aNurb.Pw(ipoint,:);
        end
        aNurb.Pw = auxPw;
        nurbs(iNurbs) = aNurb;
    end
end


%______________________________________
function [ElemNodesInfo,nonConNodes] = createElemNodesInfo(T,exteriorNodes,faceVertices)

nOfNodes = length(exteriorNodes);
nOfElements = size(T,1); 
maxNodes = max(exteriorNodes);
ElemNodesInfo = spalloc(maxNodes,2,2*nOfNodes);
num = 1:nOfNodes;
pos = ones(nOfNodes,1);
for elem = 1:nOfElements
    iFace = T(elem,end);
    faceNodes = faceVertices(iFace,:);
    nodes = T(elem,faceNodes);
    for j = [1 2]
        jnode = num(exteriorNodes == nodes(j));
        posNode = pos(jnode);
        ElemNodesInfo(nodes(j),posNode) = elem;
        pos(jnode) = posNode + 1;
    end
end
nonConNodes = exteriorNodes(pos == 2);
    
    




