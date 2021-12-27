function uint = Berkhoff1D(X1,T1,referenceElement,k1,kx,ky,z1,omega,alpha)

if ischar(alpha) && strcmpi(alpha,'rad') %Add PML to the left boundary
    withPML = true;
    
    %PML lenght = \lambda * pml_lonFactor
    pml_lonFactor = 3;
    R = 40;
    pml_n = 2;
    pml_sigma0 = R*(pml_n+1)/2;
    p = referenceElement.degree;
    h = abs(X1(p+1) - X1(1));

    %PML mesh
    L_A = pml_lonFactor*2*pi/k1(1);
    A = X1(1) - L_A;
    elems_A = round(abs(X1(1)-A)/h);
    nOfNodes_A = p*elems_A + 1;
    Tpml_A = create1Dconec(elems_A,p);
    Xpml_A = map1Dmesh(linspace(X1(1),A,nOfNodes_A)',Tpml_A,referenceElement);
    Xpml_A = flipud(Xpml_A);
    Taux = T1 + length(Xpml_A) - 1;
    Tpml_A(end) = Taux(1);

    %Absorption parameter in the PML
    sigma_A = pml_sigma0*abs(((Xpml_A-Xpml_A(end))/L_A)).^pml_n;

    %Extended mesh
    X = [Xpml_A(1:end-1) ; X1];
    T = [Tpml_A ; Taux];

    %Absorption parameter in the extended mesh
    sigma = zeros(size(X));
    sigma(unique(Tpml_A)) = sigma_A;

    %Other parameters in the extended mesh
    k = [k1(1)*ones(size(Xpml_A(2:end))) ; k1];
    z = [z1(1)*ones(size(Xpml_A(2:end))) ; z1];
    
else %No PML
    withPML = false;
    
    nOfNodes_A = 1;
    X = X1;
    T = T1;
    k = k1;
    z = z1;
    sigma = zeros(size(X));
end

%Celerities
ccg = celerities(omega,k,z);
nOfNodes = length(X);

%Matrices and rhs vector
[K,M] = matrices1D_Berkhoff(X,T,referenceElement,sigma,ccg,k,kx,omega);
Mat = K - M;
f = zeros(nOfNodes,1);

%Boudary condition on the left boundary
im = sqrt(-1);
if withPML
    %Natural BC
else
    %Partial reflecting boundary
    Mat(1,1) = Mat(1,1) - alpha*im*ccg(1)*k(1);
end

%Essential BC on the right boundary
index = nOfNodes;
value = exp(im*ky*X(index));
column = Mat(:,index);
column(index) = 0;
Mat(index,:) = 0;
Mat(:,index) = 0;
Mat(index,index) = 1;
f(index) = value;
f = f - value*column;

%Solution
u = Mat \ f;

%Solution in the interior domain
uint = u(nOfNodes_A:end);


























