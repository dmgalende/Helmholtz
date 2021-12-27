function  data = createBottomDataAUX(data,file)

if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
    nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
elseif strcmp(data.computation,'NEFEM')
    nameNodal = 'X';
end

if file ~= 0
    bottomFun = str2func(file(1:end-2));
    bottom = bottomFun(data.mesh.(nameNodal));
    if isfield(data,'bottom')
        data = rmfield(data,'bottom');
    end
    data.bottom.value = bottom;
    return
end

%Check PML condition
PML = [];
nOfBoundaries = numel(data.mesh.boundaryNames);
totalPML = false(1,nOfBoundaries);
boundaries = 1:nOfBoundaries;
for iboundary = 1:nOfBoundaries
    if strcmp(data.PML{5,iboundary},'on')
        PML = [PML iboundary];
    end
end
if isempty(PML)
    return
end

totalPML = boundaries(PML);
nOfPMLs = length(totalPML);
computeSmoothBottom %run smooth process and create the bottom array newBottom
data.bottom.value = newBottom;

%_________________
    function computeSmoothBottom
        
        %Mesh info
        if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
            nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};
            conec = data.mesh.(nameCon);
        elseif strcmp(data.computation,'NEFEM')
            nOfBoundaries = numel(data.mesh.boundaryNames);
            conec = data.mesh.T;
            for iboundaryNEFEM = 1:nOfBoundaries
                inameCon = data.mesh.fieldNames{data.mesh.boundaryIndex(iboundaryNEFEM)};
                conec = [conec ; data.mesh.(inameCon)];
            end
        end
        X = data.mesh.(nameNodal);
        
        %PML info
        PMLinfo = cell(1,4,nOfPMLs); %Coordinates of each cartesian PML boundary
        PMLlon = zeros(nOfPMLs,1); %Length of each PML area
        for i = 1:nOfPMLs
            iPML = data.PML{2,totalPML(i)};
            PMLinfo(:,:,i) = iPML;
            lon = data.BC.parameters{totalPML(i)}{4}(3);
            PMLlon(i) = lon;
        end
        
        %----------------------------------------
        % CONSTANT BOTTOM ALONG THE PML DIRECTION
        %----------------------------------------
       
        %Nodal Connectivity
        Nlin = getNodalConnectivity(conec(:,1:3));
        
        %Initialize new bottom
        newBottom = data.bottom.value;
        
        %Loops in PMLs, PML type and PML fragment
        coords_aux = [1 2 1 2];
        PMLcoefficients_aux = [1 -1 -1 1];
        PMLpermutation_aux = [1 2 3 4 ; 4 3 1 2 ; 1 2 4 3 ; 3 4 1 2];
        tol_aux = 1e-5;
        numvec = 1:size(X,1);
        for i = 1:nOfPMLs
            for j = 1:4
                if ~isempty(PMLinfo{:,j,i})
                    
                    %Coordinates for the ith-PML of type j
                    jx = coords_aux(j);
                    jy = 3 - jx;
                    jcoef = PMLcoefficients_aux(j);
                    jperm = PMLpermutation_aux(j,:);
                    
                    for k = 1:size(PMLinfo{:,j,i},1)
                        
                        %Box for the kth-fragment of the ith-PML of type j
                        kl1 = PMLinfo{:,j,i}(k,1) - tol_aux;
                        kl2 = PMLinfo{:,j,i}(k,2) + tol_aux;
                        kl3 = PMLinfo{:,j,i}(k,3) - jcoef*tol_aux;
                        kl4 = PMLinfo{:,j,i}(k,3) + jcoef*PMLlon(i) + jcoef*tol_aux;
                        klimits = [kl1 kl2 kl3 kl4];
                        
                        %Mesh elements inclouded in this PML
                        kElemPosPML = getElemPosFromSquareBox(X,conec,klimits(jperm));
                        kT_PML = conec(kElemPosPML,:);
                        kT_PMLlin = kT_PML(:,1:3);
                        
                        %PML nodes
                        knodesPML = unique(kT_PML);
                        kcond = X(knodesPML,jy) < PMLinfo{:,j,i}(k,3) + jcoef*PMLlon(i) + tol_aux & ...
                            X(knodesPML,jy) > PMLinfo{:,j,i}(k,3) + jcoef*PMLlon(i) - tol_aux;
                        knodesinPML = knodesPML(~kcond);
                        
                        %Linear interface nodes of this PML
                        auxvec = false(size(X,1),1);
                        auxvec2 = auxvec;
                        knodesPMLlin = unique(kT_PMLlin);
                        auxvec(knodesPML(kcond)) = true;
                        auxvec2(knodesPMLlin) = true;
                        knodesiPMLlin = numvec(auxvec & auxvec2);
                        
                        %Update new bottom for this PML nodes
                        newBottom = computeNewBottom(X,conec,Nlin,knodesiPMLlin,knodesinPML,newBottom,jx);
                    end
                end
            end
        end
        
    end
%_________________


%_________________
    function bottom = computeNewBottom(X,T,N,basepoints,interpoints,bottom,jx)
        
        [xi,opos] = sort(X(basepoints,jx));
        obasepoints = basepoints(opos);
        xi_t = xi';
        x = X(interpoints,jx);
        difmat = sign(x(:,ones(1,size(xi_t,2))) - xi_t(ones(size(x,1),1),:));
        posmax = size(difmat,2);
        
        for i = 1:numel(interpoints)
            
            %Find vertex nodes
            ipos = 1; while difmat(i,ipos) >= 0 && ipos < posmax, ipos = ipos + 1; end
            v1 = obasepoints(ipos-1);
            v2 = obasepoints(ipos);
            v = [v1 v2];
            
            %Find element
            elem1 = N(v1,logical(N(v1,:)));
            elem2 = N(v2,logical(N(v2,:)));
            cond = false; j = 1; while ~any(cond), cond = elem1(j) == elem2; j = j+1; end
            elem = elem2(cond);
            
            %Find face
            fn = data.mesh.referenceElement.faceNodes(:,[1 end]);
            j = 1; while any(T(elem,fn(j,:)) ~= v) && any(T(elem,fn(j,:)) ~= fliplr(v)); j = j+1; end
            
            %Interpolate bottom
            inodes = T(elem,data.mesh.referenceElement.faceNodes(j,:));
            ix = X(inodes,jx);
            ibottom = bottom(inodes);
            p = vander(ix) \ ibottom;
            bottom(interpoints(i)) = polyval(p,x(i));
        end
        
    end
%_________________


%_________________
    function N = getNodalConnectivity(T_delaunay)
        
        nNodes = max(max(T_delaunay));
        nNodesElem = size(T_delaunay,2);
        N = zeros(nNodes,10);
        nn = ones(nNodes,1);
        for ielem = 1:size(T_delaunay,1)
            Te = T_delaunay(ielem,:);
            nn_Te = nn(Te);
            for kk = 1:nNodesElem
                N(Te(kk),nn_Te(kk)) = ielem;
            end
            nn(Te) = nn(Te) + 1;
        end
        N(:,max(nn):end) = [];
        
    end
%_________________


%_________________
    function res = getElemPosFromSquareBox(X,T,limits)
        
        X1 = X(:,1);
        X2 = X(:,2);
        
        x = X1(T);
        y = X2(T);
        
        auxmatrix = x > limits(1) & x < limits(2) & y > limits(3) & y < limits(4);
        
        res = all(auxmatrix,2);
        
    end
%_________________

end



