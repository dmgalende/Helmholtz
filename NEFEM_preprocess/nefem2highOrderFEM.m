home, clear all

coords  = load('coords_readed.mat');
nefem   = load('BellottiHarbor_NEFEM1[p4][h150].mat');

nefemfields = fieldnames(nefem);
boundaries = {};
cont = 1;
for i = 1:length(nefemfields)
    ifield = nefemfields{i};
    if length(ifield) >= 3 && strcmp(ifield(1:3),'Tb_')
        boundaries{cont} = ifield;
        cont = cont + 1;
    end
end

refElem = createReferenceElement(1,coords.elemInfo.nOfNodes);
tol = 1e-5;
for c_elem = 1:size(coords.T,1)
    for ibdry = 1:length(boundaries)
        elem = 0;
        name = boundaries{ibdry};
        for i = 1:size(nefem.(name),1)
            node = 1;
            nodecond = true;
            while nodecond && node <= 3
                nodex = nefem.X(nefem.(name)(i,node),1) >= coords.X(coords.T(c_elem,1:3),1) - tol & ...
                    nefem.X(nefem.(name)(i,node),1) <= coords.X(coords.T(c_elem,1:3),1) + tol;
                nodey = nefem.X(nefem.(name)(i,node),2) >= coords.X(coords.T(c_elem,1:3),2) - tol & ...
                    nefem.X(nefem.(name)(i,node),2) <= coords.X(coords.T(c_elem,1:3),2) + tol;
                nodecond = nodecond & any(nodex & nodey);
                node = node + 1;
            end
            if nodecond, elem = i; break, end
        end
        if elem == 0 && ibdry < length(boundaries)
            continue
        elseif elem == 0 && ibdry == length(boundaries)
            warning('Did not find a matching element in the nefem mesh')
            continue
        else
            disp(['Element ' num2str(i) ' found'])
        end

        face = 1; %nefem element has the curved boundary in the 1st face!
        Xfacetri = nefem.X(nefem.(name)(elem,coords.elemInfo.faceNodes(face,[1,end])),:);
        facebool = [0 0 0]';
        for j = 1:2
            facenodex = Xfacetri(j,1) >= coords.X(coords.T(c_elem,1:3),1) - tol &...
                Xfacetri(j,1) <= coords.X(coords.T(c_elem,1:3),1) + tol;
            facenodey = Xfacetri(j,2) >= coords.X(coords.T(c_elem,1:3),2) - tol &...
                Xfacetri(j,2) <= coords.X(coords.T(c_elem,1:3),2) + tol;
            facebool = facebool | (facenodex & facenodey);
        end
        if facebool == [1;1;0]
            face_orig = 1;
        elseif facebool == [1;0;1]
            face_orig = 3;
        elseif facebool == [0;1;1]
            face_orig = 2;
        else
            error('Did not find any matching face in the nefem mesh')
        end

        nodesface_orig = refElem.faceNodes(face_orig,:);
        nodesface = refElem.faceNodes(face,:);
        nefem.X(nefem.(name)(elem,nodesface),:) = coords.X(coords.T(c_elem,nodesface_orig),:);
        nefem.X(nefem.(name)(elem,refElem.innerNodes),:) = coords.X(coords.T(c_elem,refElem.innerNodes),:);
        
        break
    end
end

conec = [];
ini = 1;
fin = 0;
for i = 1:length(boundaries)
    iname = boundaries{i};
    ivar = zeros(size(nefem.(iname),1),size(coords.elemInfo.faceNodes,2));
    for j = 1:size(ivar,1)
        jface = 1; %nefem element has the curved boundary in the 1st face!
        ivar(j,:) = nefem.(iname)(j,coords.elemInfo.faceNodes(jface,:));
    end
    conec = [conec ; nefem.(iname)];
    nefem.(iname) = ivar;
    
    fin = fin + size(ivar,1);
    numelems = (ini:fin)';
    ini = fin + 1;
    nefem.elementFaceInfo.(iname(4:end)) = [numelems,ones(size(numelems))]; %nefem element has the curved boundary in the 1st face!
end
nefem.T = [conec ; nefem.T];

nefem = rmfield(nefem,{'nurbs' 'trimmedInfo'});


