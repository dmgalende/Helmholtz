function h = bottomFunctionShoaling_var(X)

x = X(:,1);
h = -x/100 + 10; %not constant
% h(h>10) = 10;
h(h<2) = 2;

