function [beta, beta2, dbeta, dbeta2, ddbeta, ddbeta2] = eval_KDFE_coefs(m, k, r, H0_, H1_)

if nargin == 3
    H0 = besselh(0, k*r);
    H1 = besselh(1, k*r);
else
    H0 = H0_;
    H1 = H1_;
end

dH0 = -k * H1;
dH1 = k * H0 - (1/r) * H1;
d2H0 = -k * dH1;
d2H1 = k * dH0 + (1/r^2) * H1 - (1/r) * dH1;

aux = (k*r)^m;
beta = H0 / aux;
beta2 = H1 / aux;

aux = 1 / (k*r)^(m+1);
coef = k*r*dH0 - k*m*H0;
coef2 = k*r*dH1 - k*m*H1;
dbeta = aux * coef;
dbeta2 = aux * coef2;

aux2 = -k*(m+1) / (k*r)^(m+2);
ddbeta = aux2 * coef + aux * (k*r*d2H0 + k*dH0*(1-m));
ddbeta2 = aux2 * coef2 + aux * (k*r*d2H1 + k*dH1*(1-m));
