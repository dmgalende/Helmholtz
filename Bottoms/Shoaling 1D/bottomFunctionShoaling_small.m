function h = bottomFunctionShoaling_small(X)

x = X(:,1);
h =  10 - 8* x/max(x) ; %not constant
h(h>10) = 10;
h(h<2) = 2;

