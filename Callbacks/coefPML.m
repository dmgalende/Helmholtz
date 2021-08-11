function [sigma,finalSigmaPos,stretching] = coefPML(Xmesh,nodes,boundary,n,R,lon)

X = Xmesh(nodes,:);
sigma0 = R*(n+1)/2; %sigma0 = ((n+1)/(2*lon))*R; && modified 14/05/2010
Fcoef = sqrt(-1)*sigma0/((n+1)*(lon^n));
nOfNodes = size(X,1);
coord2check = [1 2 1 2];
coef = [1 -1 -1 1];
tol = 1e-3;

sigma = zeros(nOfNodes,2);
stretching = zeros(size(X));
nodesNum = 1:nOfNodes;
sigmaPos = false(nOfNodes,2);
finalSigmaPos = cell(1,2); %Initialize the PML nodes position

for i = 1:nOfNodes
    ipos = X(i,:);
    for j = 1:4
        jtype = boundary{j};
        jcoef = coef(j);
        jcoord = coord2check(j);
        jsave = 3 - jcoord;
        jlon = jcoef*lon;
        for k = 1:size(jtype,1)
            kboundary = jtype(k,:);
            iniPML = kboundary(3);
            endPML = iniPML + jlon;
            if tol+ipos(jcoord) >= kboundary(1) && ipos(jcoord) <= kboundary(2)+tol &&...
                    tol+jcoef*ipos(jsave) >= jcoef*kboundary(3) &&...
                    jcoef*ipos(jsave) <= jcoef*endPML+tol
                
                %Polinomyal absorption parameter
                sigma(i,jsave) = sigma0*abs(((ipos(jsave)-endPML)/lon))^n;
                
                %Coordinate stretching
                stretching(i,jsave) = Fcoef*(abs(ipos(jsave)-endPML)^(n+1));
                
                %Marked node as PML one
                sigmaPos(i,jsave) = true;
                break
            end
        end
    end
end

%PML nodes numbering
for i = [1 2]
    finalSigmaPos{i} = nodesNum(sigmaPos(:,i));
end
