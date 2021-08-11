function h = bottomFunctionRefraction(X)

x = X(:,1);
h = -0.05*x + 35; %not constant
pos = find(h > 35);
h(pos) = 35;

