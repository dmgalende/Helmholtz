function Xm = map1Dmesh(X,T,referenceElement)

%Mapping a 1D high-order uniform mesh {X,T} into a non-uniform mesh (fekete type) Xm

p = referenceElement.degree;
Xm = X;
if p > 2
    refcoord = referenceElement.NodesCoord1d;
    shapeFunctions = (1/2)*[1-refcoord,1+refcoord];
    for ielem = 1:size(T,1)
        Te = T(ielem,:);
        Xlin = X(Te([1,end]));
        Xm(Te) = shapeFunctions*Xlin;
    end
end