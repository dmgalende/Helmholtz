function varargout = plotSolution_NEFEM(mesh,u)

%Boundary info
meshFields = fieldnames(mesh);
T = [];
for i = 1:length(meshFields)
    name = meshFields{i};
    if length(name) > 3 && strcmpi(name(1:3),'Tb_'), 
        T = [T ; mesh.(name)];
    end
end

%Plot solution with FEM functions
h = plotSolution(mesh.X,[T ; mesh.T],u,mesh.referenceElement);

%Output variable
if ~nargout
    varargout = [];
elseif nargout == 1
    varargout = {h};
end
