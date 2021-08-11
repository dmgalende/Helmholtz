function [ccg,c,cg] = celerities(omega,k,h)

g = 9.81;
c = omega./k;
cg = (g/(2*omega))*(tanh(k.*h) + k.*h.*(sech(k.*h)).^2);
ccg = c.*cg;
