function h = Example2(X)

pos = X(:,2) <= 2;
y = X(pos,2);
maxValue = max(y);
minValue = min(y);
hmax = 0.3;
hmin = 0.01;
h = hmax*ones(size(X,1),1);
h(pos) = ((y - minValue)/(maxValue - minValue))*(hmax - hmin) + hmin;