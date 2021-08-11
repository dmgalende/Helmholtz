function coordRef = nefemInterp2DAdaptedNodesElement(coordRefStr, vertCoord, aNurbs, u1, u2)
%
%coordRef = nefemInterp2DAdaptedNodesElement(coordRefStr, vertCoord,
%                                       aNurbs, u1, u2)
%
% Compute nodal distribution adpated to the NURBS side
% Same as nefemQuad2DElementLocalCoordinates avoiding weights computation
%

TOL = 1e-8;

vertexParam = [0 1];
nOfNodes = size(coordRefStr, 1);
coordRef = zeros(nOfNodes,2);
for iNode = 1:nOfNodes
    if coordRefStr(iNode,2) == 1
        point_rst = [-1,1];
    else
        point_rst = mapXiEtaZeta2rst(coordRefStr(iNode,:));
    end
    % On [0,1]^2
    point_rst = (point_rst+1)/2;

    u = point_rst(1)*(u2 - u1) + u1;
    if abs(u) < TOL
        u = 0;
    elseif abs(u-u1) < TOL
        u = u1;
    elseif abs(u-u2) < TOL
        u = u2;   
    end
    v = point_rst(2);

    pt = nurbsCurvePoint(aNurbs,u);
    phi01 = inverseLinearMapping(vertCoord,pt(1:2));
    phi = (phi01+1)/2;
    
    coordRef(iNode,:) = (1-v)*phi  + v*vertexParam;
end

% On [-1 -1; 1 -1; -1 1]
coordRef = 2*coordRef - 1;

