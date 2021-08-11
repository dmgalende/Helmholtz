function h = ellipticShoal(X)

x0 = 0;
y0 = 0;
angle = 20;
x = X(:,1);
y = X(:,2);
xp = (x-x0)*cos(pi/180*angle) + (y-y0)*sin(pi/180*angle);
yp = (y-y0)*cos(pi/180*angle) + (x-x0)*sin(pi/180*angle);
h = 0.45*ones(size(X,1),1);
% bottom on the linear shoal 
h(yp<5.84) = 0.45-0.02*(5.84-yp(yp<5.84));
% bottom on the elliptic shoal
ellipse = (xp/4).^2+(yp/3).^2;
h(ellipse<1) = h(ellipse<1) +0.3-0.5*sqrt(1-(xp(ellipse<1)/5).^2-(yp(ellipse<1)/3.75).^2);
