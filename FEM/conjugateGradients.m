function [x,flag,k,errork] = conjugateGradients(A,b,x0,nmax,tol)

tol_default = 1e-6;
nmax_default = 100;
tolmu = 1e-15;
switch nargin
    case 3
        tol = tol_default;
        nmax = nmax_default;
    case 4
        tol = tol_default;
end
nOfComputesRes = 0.3*nmax;

r0 = A*x0 - b;
normr0 = norm(r0);
p0 = r0;
rho0 = sdot(r0,r0);

errork = zeros(nmax,1);
x = x0;
convergence = false;
k = 0;
rkvalues = ceil(linspace(k,nmax,nOfComputesRes));
rkvalue = rkvalues(2);
i = 3;
while ~convergence && k < nmax
    
    qk = A*p0;
    muk = sdot(p0,qk);
%     if abs(muk) < tolmu, flag = 2; return, end
    alphak = -rho0/muk;
    x = x0 + alphak*p0;
    rk = r0 + alphak*qk;
%     rk = A*x - b;
%     if k == rkvalue, rk = A*x - b; rkvalue = rkvalues(i); i = i + 1; else rk = r0 + alphak*qk; end
    rhok = sdot(rk,rk);
    betak = rhok/rho0;
    pk = rk + betak*p0;
    
    k = k + 1;
    
    errork(k) = norm(rk)/normr0;
    convergence = errork(k) <= tol;
    
    x0 = x;
    r0 = rk;
    p0 = pk;
    rho0 = rhok;
    
end

if convergence, flag = 1; errork = errork(1:k); else flag = 0; end


function res = sdot(a,b)

% res = a'*b;
res = a.'*b;


    