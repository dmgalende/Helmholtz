function [u0,u0gradStruct,der2DShapeFunOn1D,Text,Tbext,Tint,inPMLelems,elementFaceInfo_boundaryPML] =...
    computeIncidentWaveField(data,kvec,omega,handle)

%% INCIDENT WAVE PREPROCESS

%Mesh & parameters
X = data.mesh.X;
T = data.mesh.T;
ne = size(T,1);
nOfNodes = size(X,1);
referenceElement = data.mesh.referenceElement;
p = referenceElement.degree;
kx = kvec(1);
ky = kvec(2);
im = sqrt(-1);

%Get PML and problem type
pmlflag = false;
isScat = false; %Scattering problem type (no harbors)
for iboundary = 1:numel(data.mesh.boundaryNames)
    pml_aux = data.PML(:,iboundary);
    if strcmp(pml_aux{5},'on')
        pml = pml_aux;
        pml_params = data.BC.parameters{iboundary}{4};
        if ~pmlflag, pmlflag = true;
        else
            setOutput({'Up to now, routine computeIncidentWaveField.m only accepts ONE PML'},handle)
            error('error in Berkhoff GUI')
        end
        isBoundedPML = ~isempty(pml{2}{1}); %Check if PML has a down boundary
        if isBoundedPML, isScat = true; break, end %Force to be a scattering problem type for closed PMLs
    else
        if strcmpi(data.mesh.boundaryNames{iboundary},'ext')
            cond_ext = data.BC.values(iboundary);
            if any(cond_ext == [2,3]) %Ext boundary always artificial
                alpha = 'rad';
            elseif cond_ext == 1
                isScat = true; %Force to be a scattering problem type for symmetry conditions
            end
        end
    end
end

setOutput({'    Getting PML interface data...'},handle)

%PML connectivity
[out_e,~,Nconec] = createOutOfPMLConnectivity(data);
Tint = T(out_e,:);
inPMLelems = setdiff(1:ne,out_e);
Text = T(inPMLelems,:);
nodesPML = unique(Text);
nodesINT = setdiff(1:nOfNodes,nodesPML);
cont_nodes = intersect(nodesPML,unique(Tint));

%Left interface
node_aux = pml{3}{4}(1);
xL_aux = X(node_aux,1) + pml_params(3);
leftNodes = cont_nodes(X(cont_nodes,1) < xL_aux + abs(1e-5*xL_aux));
xL = X(leftNodes(1),1);
yL = min(X(leftNodes,2));
dL = data.bottom.value(leftNodes);

%Right interface
node_aux = pml{3}{2}(1);
xR_aux = X(node_aux,1) - pml_params(3);
rightNodes = cont_nodes(X(cont_nodes,1) > xR_aux - abs(1e-5*xR_aux));
xR = X(rightNodes(1),1);
yR = min(X(rightNodes,2));
dR = data.bottom.value(rightNodes);

%Constant bathymetry area
y0 = max(X(rightNodes,2));
y0_L = max(X(leftNodes,2));
if abs((y0-y0_L)/y0_L) > 1e-2
    error('Constant bathymetry area on right and left PML interfaces dont coincide')
end
restNodes = setdiff(cont_nodes,[leftNodes ; rightNodes]);

%Up interface
upNodes = restNodes(X(restNodes,2) > y0 - abs(1e-5*y0));

if ~isScat
    setOutput({'    Boundary EXT is artificial'},handle);
    
    %There are not down PML inteface
    downNodes = [];

    %1D solution in the y direction on the LEFT and RIGHT interfaces evaluated at all coordinates X(:,2)
    setOutput({'    Computing 1D problems...'},handle);
    factor_h = 1/20;
    min_p = 6;
    uL = main_Berkhoff1D(X,nodesPML,y0_L,yL,dL,leftNodes,referenceElement,...
        factor_h,min_p,omega,kx,ky,alpha);
    uR = main_Berkhoff1D(X,nodesPML,y0,yR,dR,rightNodes,referenceElement,...
        factor_h,min_p,omega,kx,ky,alpha);

    %Splines in the x direction for LEFT and RIGHT interfaces evaluated at all coordinates X(:,1)
    tauL = 0.5;
    tauR = 0.5;
    [sL,s0,sR] = interpolationSplines(tauL,tauR,X(:,1),xL,xR);

    %Incident wave field over the computational domain
    u0 = (sL.*uL + s0.*exp(im*ky*X(:,2)) + sR.*uR) .* exp(im*kx*X(:,1));

else
    
    %Down nodes
    if isBoundedPML
        downNodes = setdiff(restNodes,upNodes);
    else
        downNodes = [];
    end
    
    %Incident wave for constant bathymetry
    setOutput({'    Incident wave for constant bathymetry'},handle);
    u0 = data.ip.amplitude * exp(im*(kx*X(:,1) + ky*X(:,2)));
end

%The values of u0 outside the PML are not needed for the formulation with the total wave
if ~data.staticCondensation, u0(nodesINT) = 0; 
else u0 = data.ip.amplitude * exp(im*(kx*X(:,1) + ky*X(:,2))); 
end %STATIC COND. SHOULD BE FIXED WITH TOTAL WAVE FORM.!!! Then remove this condition and let u0(nodesINT) = 0

%% COMPUTING GRADIENTS ON THE 1D GAUSS POINTS OF THE MESH BOUNDARY

setOutput({'    Computing gradients on 1D integration points...'},handle);
u0gradStruct = struct();
der2DShapeFunOn1D = struct();
Text_nodes = [];
prev = 1;
post = 0;
for iboundary = 1:numel(data.mesh.boundaryNames)
    boundaryName = data.mesh.boundaryNames{iboundary};
    Tb = data.mesh.(['Tb_' boundaryName]);
    
    if ~isScat
        [u0gradStruct.(boundaryName),AUXder2DShapeFunOn1D] = compute2DScalarFieldGradOnOneDimIP(...
            u0,...
            X,T,...
            Tb,...
            data.mesh.elementFaceInfo.(boundaryName),...
            data.mesh.referenceElement);
    else
        [u0gradStruct.(boundaryName),AUXder2DShapeFunOn1D] = compute2DIncidentWaveFieldGradOnOneDimIP(...
            X,T,...
            Tb,...
            data.mesh.elementFaceInfo.(boundaryName),...
            data.mesh.referenceElement,...
            kvec,...
            data.ip.amplitude);
    end
    
    if any(strcmpi(boundaryName,{'ext' 'pml'}))
        inelems = size(Tb,1);
        post = post + inelems;
        Text_nodes = [Text_nodes ; Tb];
        elementFaceInfo_boundaryPML(prev:post,:) = data.mesh.elementFaceInfo.(boundaryName);
        u0gradStruct.extTb(:,:,prev:post) = u0gradStruct.(boundaryName);
        der2DShapeFunOn1D.extTb(:,:,:,prev:post) = AUXder2DShapeFunOn1D;
        prev = post + 1;
    end
end

%% PML INTERFACE DATA

%1D mesh of the PML interface
[~,oposL] = sort(X(leftNodes,2),'ascend');
[~,oposR] = sort(X(rightNodes,2),'descend');
[~,oposU] = sort(X(upNodes,1),'ascend');
[~,oposD] = sort(X(downNodes,1),'descend');
list_cont_nodes = [leftNodes(oposL) ; upNodes(oposU) ; rightNodes(oposR) ; downNodes(oposD)];
if isBoundedPML, cont_elems = (length(list_cont_nodes)-p)/p;
else cont_elems = (length(list_cont_nodes)-1)/p; 
end
Tcont_nodes = create1Dconec(cont_elems,p,list_cont_nodes);

%elementFaceInfo and gradients on the PML interface
elementFaceInfo_cont = computeElementFaceInfo(Tcont_nodes,Nconec,T,inPMLelems);
post = post + size(Tcont_nodes,1);
if ~isScat
    [u0gradStruct.extTb(:,:,prev:post),der2DShapeFunOn1D.extTb(:,:,:,prev:post)] =...
        compute2DScalarFieldGradOnOneDimIP(...
            u0,...
            X,T,...
            Tcont_nodes,...
            elementFaceInfo_cont,...
            data.mesh.referenceElement);
else
    [u0gradStruct.extTb(:,:,prev:post),der2DShapeFunOn1D.extTb(:,:,:,prev:post)] =...
        compute2DIncidentWaveFieldGradOnOneDimIP(...
            X,T,...
            Tcont_nodes,...
            elementFaceInfo_cont,...
            data.mesh.referenceElement,...
            kvec,...
            data.ip.amplitude);
end

%1D mesh of the bounded boundary of the exterior domain (PML area)
elementFaceInfo_boundaryPML(prev:post,:) = elementFaceInfo_cont;
Tbext = [Text_nodes ; Tcont_nodes]; %Tcont_nodes in last order!




