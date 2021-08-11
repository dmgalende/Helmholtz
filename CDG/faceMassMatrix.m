function [M_face] = faceMassMatrix(theReferenceElement)

nOfFaceNodes = numel(theReferenceElement.faceNodes1d);
IPw = theReferenceElement.IPweights1d; 
N = theReferenceElement.N1d;
%Number of Gauss points
ngauss = length(IPw);

M_face = zeros(nOfFaceNodes,nOfFaceNodes);

for g = 1:ngauss
    %Shape functions at the current integration point
    N_g = N(g,:);
    %Integration weight
    dL = IPw(g);
    M_face = M_face +  N_g' * N_g *  dL;
end

