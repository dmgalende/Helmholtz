function data = Initialization(data)
%
% [U,t] = maxwellInitialization(mesh,theReferenceElement)
%
% Initilization
if any(strcmp(data.PML(5,:),'on'))
    nOfComponents = 5;
else
    nOfComponents = 4;
end
nElem = size(data.mesh.T,1);
nOfElementNodes =  size(data.mesh.T,2);
data.U = zeros(nOfElementNodes,nOfComponents,nElem);
data.t = 0;
