function h = bottomFunHomma_modRad(X)

hb = 30;
ha = 3;
a = 600;
b = 1800;

A = (hb-ha)/(b^2-a^2);
B = ha - A*a^2;

h = zeros(size(X,1),1);
for i = 1:length(X(:,1))
    [theta,r] = cart2pol(X(i,1),X(i,2));
    if r <= b
        h(i) = A*r^2 + B;
    else
        h(i) = hb;
    end
end


