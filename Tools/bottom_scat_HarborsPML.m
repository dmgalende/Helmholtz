function d = bottom_scat_HarborsPML(X)

k = 0.0576; 
g = 9.81; 
w = 2*pi/10;

d = zeros(size(X,1),1);
d(:) = (1/k)*atanh(w^2/(k*g));