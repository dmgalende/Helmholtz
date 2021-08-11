function [M_face] = faceMassMatrix_x(Nf1,Nf2,refEl_1,refEl_2)

nOfFaceNodes_1 = numel(refEl_1.faceNodes1d);
nOfFaceNodes_2 = numel(refEl_2.faceNodes1d);
IPw = refEl_1.IPweights1d; 
%Number of Gauss points
ngauss = length(IPw);

M_face = zeros(nOfFaceNodes_1,nOfFaceNodes_2);

for g = 1:ngauss
    %Shape functions at the current integration point
    Nf1_g = Nf1(g,:);
    Nf2_g = Nf2(g,:);
    %Integration weight
    dL = IPw(g);
    M_face = M_face +  Nf1_g' * Nf2_g *  dL;
end

