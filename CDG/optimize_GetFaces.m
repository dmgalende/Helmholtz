load('/media/My Passport/Berkhoff_GUI/Saves/ScatteringCircle/Comparison_CDG-FEM/data__T_0.2_kh_0.75_PML_0.2_P1_doubled.mat')

[intFaces,extFaces] = GetFaces(data.mesh.T(:,1:3));