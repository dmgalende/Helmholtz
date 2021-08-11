function edgeNodesOrdered = postprocessComputeEdgeNodesOrdered(p)


nNodesEdge = p + 1;
edgeNodesOrdered = 1:nNodesEdge;
for i = 2:p+1
    edgeNodesOrdered = [edgeNodesOrdered, edgeNodesOrdered(end)+nNodesEdge-i+1];
end
for i = 2:p+1
    edgeNodesOrdered = [edgeNodesOrdered, edgeNodesOrdered(end)-i];
end
