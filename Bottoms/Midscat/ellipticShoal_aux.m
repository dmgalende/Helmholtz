function h = ellipticShoal_aux(X)

Rc = 1; %radius of the circle
a = 0.154305777009007;
b = 0.115551367937199;
angle = 20;

x = X(:,1) / Rc; y = X(:,2) / Rc;
tol = 1e-5;
pos1 = x <= 2+tol & x >= -2-tol;
pos2 = x < -2-tol;
pos3 = ~(pos1 | pos2);

h = zeros(size(X,1),1);
h(pos1) = a*(y(pos1)*cos(pi/180*angle) + x(pos1)*sin(pi/180*angle)) + b;
h(pos2) = a*(y(pos2)*cos(pi/180*angle) - 2*sin(pi/180*angle)) + b;
h(pos3) = a*(y(pos3)*cos(pi/180*angle) + 2*sin(pi/180*angle)) + b;
h(h>0.3) = 0.3;
h(h<0.01) = 0.01;
