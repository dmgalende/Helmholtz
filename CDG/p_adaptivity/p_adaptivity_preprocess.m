function data = p_adaptivity_preprocess(data)
omega = 2*pi/data.ip.period;
elem_p = data.p_adaptivity.elem_p;
nOfElements = length(elem_p);
p_base = data.mesh.referenceElement.degree;
% create check for all the elements:
% 1 refined elements; 0 base elements
% true if p_elem ~= p_base | false if p_elem == p_base
elem_check = elem_p ~= p_base;

% create the reference elements for the refined elements
check_p = false(20,1);
check_p(p_base) = true;
elements = 1:nOfElements;
nOfRefNodes = 0.5*sum((elem_p(elem_check)+1).*(elem_p(elem_check)+2));
coord = zeros(nOfRefNodes,2);
h = zeros(nOfRefNodes,1);
aux_ind = 0;
for ielem = elements(elem_check)
    p = elem_p(ielem);
    nOfElementNodes = (p+1)*(p+2)/2;
    if ~check_p(p)
        refEl = createReferenceElement(1,nOfElementNodes);
        data.p_adaptivity.referenceElements.(['P' num2str(p)]) = refEl;

        % shape functions in the faces gauss points
        Nx_face = zeros(numel(refEl.IPcoordinates1d),nOfElementNodes,3);
        Ny_face = Nx_face ;
        for iFace = 1:3
            coordGaussPointsFaceRefEl = one2twoDmapping(iFace, refEl);
            [shapeFun, shapeFunDer] = evalContinuousShapeFunctionsAtTriPoints(coordGaussPointsFaceRefEl,...
                refEl.degree, refEl.NodesCoord);
            Nx_face(:,:,iFace) = shapeFunDer(:,:,1)';
            Ny_face(:,:,iFace) = shapeFunDer(:,:,2)';
        end
        fSP.(['P' num2str(p)]).Nx_face = Nx_face;
        fSP.(['P' num2str(p)]).Ny_face = Ny_face;

        % check p
        check_p(p) = true;
        disp(['P' num2str(p) ' reference element created'])
    end

    % coordinates - ccg - k  in the refined elements
    aux_ind = aux_ind(end) + (1:nOfElementNodes);
    Te_lin = data.mesh.T(ielem,:);
    refEl = data.p_adaptivity.referenceElements.(['P' num2str(p)]);
    local_coord = refEl.NodesCoord;
    Xe = linearMapping(data.mesh.X(Te_lin,:),local_coord);
    coord(aux_ind,:) = Xe;
    h_vertices = data.bottom.value(Te_lin);
    h(aux_ind) = planarInterpolation(Xe(1:3,1),Xe(1:3,2),h_vertices,Xe(:,1),Xe(:,2));
end

% different p
different_p = find(check_p);

% ccg and wavenumber for refined elements
waveNumber = computeWaveNumber(omega,h);
ccg = celerities(omega,waveNumber,h);

% pml parameters in the refined elements
aux_boundaries = 1:numel(data.mesh.boundaryNames);
pml_index = aux_boundaries(strcmp(data.PML(5,:),'on'));
if ~isempty(pml_index)
    sigma = coefPML(coord,data.PML{2,pml_index},data.BC.parameters{pml_index}{4}(1),...
        data.BC.parameters{pml_index}{4}(2),data.BC.parameters{pml_index}{4}(3));
    data.p_adaptivity.sigma = sigma;
else
    data.p_adaptivity.sigma = zeros(nOfRefNodes,2);
end

% creating shape function structure for face loop
for ip = 1:numel(different_p)
    p = different_p(ip);
    if ip==1
        refEl = data.mesh.referenceElement;
    else
        refEl = data.p_adaptivity.referenceElements.(['P' num2str(p)]);
    end
    
    % face mass matrix
    fSP.(['P' num2str(p)]).faceMassMatrix = faceMassMatrix(refEl);
    
    p_sup = different_p(different_p>p);
    for ips = 1:numel(p_sup)
        ps = p_sup(ips);
        % higher order reference element
        refEl_s = data.p_adaptivity.referenceElements.(['P' num2str(ps)]);
        Nx_face = zeros(numel(refEl_s.IPcoordinates1d),size(refEl.NodesCoord,1),3);
        Ny_face = Nx_face ;
        for iFace = 1:3
        % faces gauss points for the element of higher order 
        coordGaussPointsFaceRefEl = one2twoDmapping(iFace, refEl_s); 
        [shapeFun, shapeFunDer] = evalContinuousShapeFunctionsAtTriPoints(coordGaussPointsFaceRefEl,...
            refEl.degree, refEl.NodesCoord);
        N_face = shapeFun'; 
        N_face = N_face(:,refEl.faceNodes(iFace,:));
        Nx_face(:,:,iFace) = shapeFunDer(:,:,1)';
        Ny_face(:,:,iFace) = shapeFunDer(:,:,2)';
        end
        fSP.(['P' num2str(ps)]).(['faceMassMatrix_x_p' num2str(p)] ) =...
            faceMassMatrix_x(refEl_s.N1d,N_face,refEl_s,refEl);
        fSP.(['P' num2str(p)]).(['N_face_x_p' num2str(ps)]) = N_face;
        fSP.(['P' num2str(p)]).(['Nx_face_x_p' num2str(ps)]) = Nx_face;
        fSP.(['P' num2str(p)]).(['Ny_face_x_p' num2str(ps)]) = Ny_face;
    end
end
data.p_adaptivity.fSP = fSP;
data.p_adaptivity.different_p = different_p;
data.p_adaptivity.elem_check = elem_check;
data.p_adaptivity.coord = coord;
data.p_adaptivity.waveNumber = waveNumber;
data.p_adaptivity.ccg = ccg;
