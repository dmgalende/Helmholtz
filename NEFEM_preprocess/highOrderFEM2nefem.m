function mesh_nefem = highOrderFEM2nefem(mesh)

mesh_fem = load(mesh);
mesh_nefem = mesh_fem;

fields = fieldnames(mesh_fem);
for i = 1:length(fields)
    ifield = fields{i};
    if length(ifield) > 3 && strcmpi(ifield(1:3),'Tb_')
        mesh_nefem.(ifield) = mesh_fem.elementFaceInfo.(ifield(4:end));
    end
end

mesh_nefem = rmfield(mesh_nefem,'elementFaceInfo');