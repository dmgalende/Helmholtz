function uint = main_Berkhoff1D(X_2D,y0,yc,bottom_interface,nodes_interface,...
    reference_2D,factor_h,min_p,omega,kx,ky,alpha)

%% MESH AND PARAMETERS FOR y in [yc,y0] 

%Mesh size & reference element
nOfNodes_2D = size(X_2D,1);
p_2D = reference_2D.degree;
[oy_interface,opos_interface] = sort(X_2D(nodes_interface,2),'ascend');
onodes_interface = nodes_interface(opos_interface);
h = abs(oy_interface(1) - oy_interface(p_2D+1)) * factor_h;
p = max(min_p,p_2D);
referenceElement = createReferenceElement(1,(p+1)*(p+2)/2);

%Mesh
nElems = round((y0-yc)/h);
nNodes = p*nElems + 1;
T = create1Dconec(nElems,p);
X = map1Dmesh(linspace(yc,y0,nNodes)',T,referenceElement);

%Interpolated bottom
obottom_interface = bottom_interface(opos_interface);
oX_2D_interface = X_2D(onodes_interface,2);
nNodes_2D_interface = size(oX_2D_interface,1);
T_2D_interface = create1Dconec((nNodes_2D_interface-1)/p_2D,p_2D);
z = interpolateSolution1D(obottom_interface,oX_2D_interface,T_2D_interface,reference_2D,X);

%Wave number
k = computeWaveNumber(omega,z);

%% DISCRETIZATION AND SOLUTION

%Solve Berkhoff equation for y in (ym,y0)
u = Berkhoff1D(X,T,referenceElement,k,kx,ky,z,omega,alpha);

%Ordered positions for coordinate y of the 2D mesh
[oy_2D,opos_2D] = sort(X_2D(:,2),'ascend');
opos_oy_extension = oy_2D >= y0;
opos_oy_interpolate = oy_2D >= yc & ~opos_oy_extension;
opos_oy_constant = ~(opos_oy_extension | opos_oy_interpolate);

%Interpolate solution on those yc =< y <= y0 of the 2D mesh
uint = zeros(nOfNodes_2D,1);
oy_interpolate = oy_2D(opos_oy_interpolate);
uint(opos_2D(opos_oy_interpolate)) = interpolateSolution1D(u,X,T,referenceElement,oy_interpolate);

%Prescribed value for those y < yc of the 2D mesh
uint(opos_2D(opos_oy_constant)) = u(1);

%Prescribed value for those y >= y0 of the 2D mesh
oy_extension = oy_2D(opos_oy_extension);
uint(opos_2D(opos_oy_extension)) = exp(sqrt(-1)*ky*oy_extension);







