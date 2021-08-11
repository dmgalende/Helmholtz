function k = computeWaveNumber(omega,bottom)

%Initial guess
g = 9.81;
if ~isscalar(omega), k0 = omega(2)*ones(size(bottom)); omega = omega(1); else k0 = omega./sqrt(g*bottom); end

%Stop parameters
nmax = 100;
tol = 1e-5;

%Pure N-R algorithm
convergence = false;
i = 0;
while ~convergence && i < nmax
    
    f = omega^2 - g*k0.*tanh(k0.*bottom);
    df = -g*tanh(k0.*bottom) - g*k0.*(1 - (tanh(k0.*bottom).^2)).*bottom;
    k = k0 - f./df;
    
    i = i + 1;
    
    errork = max(abs((k - k0)./k));
    convergence = errork <= tol;
    
    k0 = k;
    
end
k = abs(k);

if ~convergence
    warning('Newton-Raphson did not converge!!!')
end