function h = ellipticShoal(X)

x0 = 10;
y0 = 11;
angle = 180;%180-20;

x = X(:,1);
y = X(:,2);

xp = (y-y0)*sin(pi/180*angle) + (x-x0)*cos(pi/180*angle);
yp = (y-y0)*cos(pi/180*angle) - (x-x0)*sin(pi/180*angle);

h = 0.45*ones(size(X,1),1);

% bottom on the linear shoal
h(yp > -5.82) = 0.45 - 0.02*(5.82 + yp(yp > -5.82));

% bottom on the elliptic shoal
ellipse = (xp/4).^2 + (yp/3).^2;
h(ellipse < 1) = h(ellipse < 1) + 0.3 - 0.5*...
                 sqrt(1 - (xp(ellipse<1)/5).^2 - (yp(ellipse<1)/3.75).^2);
             
%Correction max depth
h(h < 0.1 & ellipse > 1) = 0.1;

