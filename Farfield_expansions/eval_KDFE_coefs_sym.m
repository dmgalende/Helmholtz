function [beta_e, beta2_e, dbeta_e, dbeta2_e, ddbeta_e, ddbeta2_e] = eval_KDFE_coefs_sym(m, k, R, H0_, H1_)

r = sym('r');
H0 = besselh(0, k*r);
H1 = besselh(1, k*r);

aux = (k*r)^m;
beta = H0 / aux;
beta2 = H1 / aux;

dbeta = diff(beta);
dbeta2 = diff(beta2);

ddbeta = diff(dbeta);
ddbeta2 = diff(dbeta2);

r = R;
beta_e = eval(beta);
beta2_e = eval(beta2);
dbeta_e = eval(dbeta);
dbeta2_e = eval(dbeta2);
ddbeta_e = eval(ddbeta);
ddbeta2_e = eval(ddbeta2);