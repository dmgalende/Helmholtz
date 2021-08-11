function h = Example3(X)

x = X(:,1);
maxValue = 1800;
minValue = -1800;
hmax = 5;
hmin = 30;
h = ((x - minValue)/(maxValue - minValue))*(hmax - hmin) + hmin;
