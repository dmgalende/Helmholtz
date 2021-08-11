function h = breakwater(X)

epsilon = 1/20;
hmin = 0.0015;
hmax = 0.305;

h = epsilon*X(:,2) + hmin;
h(h > hmax) = hmax;
