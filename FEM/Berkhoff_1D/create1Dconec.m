function T = create1Dconec(nOfElems,p,list_nodes)

T = zeros(nOfElems,p+1);
cont = 1;
if nargin == 2
    for i = 1:nOfElems
        T(i,:) = cont:cont+p;
        cont = cont+p;
    end
else
    for i = 1:nOfElems
        T(i,:) = list_nodes(cont:cont+p);
        cont = cont+p;
    end
end