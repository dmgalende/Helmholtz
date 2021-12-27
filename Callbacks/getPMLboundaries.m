function [PMLposNodes,Telem] = getPMLboundaries(X,T,referenceElement,handle,nurbs,trimmedInfo)

nOfElements = size(T,1);

elem1 = [];
elem2 = [];
elem3 = [];
elem4 = [];

tolrad = 0.5*pi/180; %0.5 degree of tolerance to detect a x or y boundary

if nargin == 4 %FEM & CDG
    Nxi = referenceElement.N1dxi(1,:);
    for ielem = 1:nOfElements
        Xe = X(T(ielem,:),:);
        xyDer = Nxi*Xe;
        xyDerNorm = norm(xyDer);
        t = xyDer/xyDerNorm;
        n = [t(2) -t(1)];
        theta = acos([1 0]*n');
        if theta <= tolrad
            elem2 = [elem2 ielem];
        elseif theta >= pi-tolrad
            elem4 = [elem4 ielem];
        elseif theta >= pi/2-tolrad && theta <= pi/2+tolrad && n(2) < 0
            elem1 = [elem1 ielem];
        elseif theta >= pi/2-tolrad && theta <= pi/2+tolrad && n(2) > 0
            elem3 = [elem3 ielem];
        else
            msje = {'??? The PML boundary doesnt seem to be a cartesian one'};
            setOutput(msje,handle), error('Error in Berkhoff GUI')
        end
    end
    
else %NEFEM
    for ielem = 1:nOfElements
        u1 = trimmedInfo(ielem).trim(1);
        u2 = trimmedInfo(ielem).trim(2);
        idNurbs = trimmedInfo(ielem).idNurbs;
        aNurbs = nurbs(idNurbs);
        igauss = nefemQuad2DFaceParametricSpace(aNurbs,u1,u2,1);
        Cu = nurbsCurveDerivPoint(aNurbs,igauss);
        jacobianNurbs = norm(Cu(1:2));
        n = [Cu(2) -Cu(1)]/jacobianNurbs;
        theta = acos([1 0]*n');
        if theta <= tolrad
            elem2 = [elem2 ielem];
        elseif theta >= pi-tolrad
            elem4 = [elem4 ielem];
        elseif theta >= pi/2-tolrad && theta <= pi/2+tolrad && n(2) < 0
            elem1 = [elem1 ielem];
        elseif theta >= pi/2-tolrad && theta <= pi/2+tolrad && n(2) > 0
            elem3 = [elem3 ielem];
        else
            msje = {'??? The PML boundary doesnt seem to be a cartesian one'};
            setOutput(msje,handle), error('Error in Berkhoff GUI')
        end
    end
    face = referenceElement.faceNodes(1,:);
    T = T(:,face);
end
Telem = {T(elem1,:) T(elem2,:) T(elem3,:) T(elem4,:)};

coord2store = [1 2 1 2];
elem = {elem1 elem2 elem3 elem4};

PMLposNodes = cell(1,4);

for i = 1:4
    Te_aux = [T(elem{i},1); T(elem{i},end)];
    Te = eliminateConectivities(Te_aux);
    Xe_aux = orderPosition(X,Te,3-coord2store(i));
    Xe = eliminateRepeatedNodes(Xe_aux); 
    PMLposNodes{i} = zeros(size(Xe,1)/2,3);
    aux = 1;
    for k = 1:size(Xe,1)/2
        PMLposNodes{i}(k,:) = [sort(Xe(aux:aux+1,coord2store(i)))'...
            Xe(aux,3-coord2store(i))];
        aux = aux+2;
    end
end


%____________________________________________________________________
function finalT = eliminateConectivities(T)

T = sort(T,'ascend');
markedRepeatedNodes = [];
nOfNodes = length(T);
for j = 1:nOfNodes
    if (j < nOfNodes && T(j) == T(j+1)) || (j > 1 && T(j) == T(j-1))
        markedRepeatedNodes = [markedRepeatedNodes j];
    end
end
T(markedRepeatedNodes) = [];
finalT = T;


%____________________________________________________________________
function finalX = eliminateRepeatedNodes(X)

markedRepeatedNodes = [];
nOfNodes = length(X);
for j = 1:nOfNodes-1
    if all(X(j,:) == X(j+1,:)) 
        markedRepeatedNodes = [markedRepeatedNodes j j+1];
    end
end
X(markedRepeatedNodes,:) = [];
finalX = X;


%____________________________________________________________________
function Xe = orderPosition(X,T,coord)

Xe_aux1 = X(T,coord);
[Xe_1,pos] = sort(Xe_aux1,'ascend');
Xe_aux2 = X(T,3-coord);
Xe = zeros(size(Xe_aux1,1),2);
Xe(:,coord) = Xe_1;
Xe_2 = Xe_aux2(pos);
breakPoint = [];
for k = 1:length(Xe_1)-1
    if Xe_1(k) ~= Xe_1(k+1)
        breakPoint = [breakPoint k];
    end
end
breakPoint = [breakPoint k+1];
aux = 1;
for u = 1:length(breakPoint)
    Xe(aux:breakPoint(u),3-coord) = sort(Xe_2(aux:breakPoint(u)),'ascend');
    aux = breakPoint(u) + 1;
end

