function ruledConnectivity = createRuledConnectivity(X,T,limits)

rule = ['x >' num2str(limits(1)) '& x <' num2str(limits(2)) ...
    '& y >' num2str(limits(3)) '& y <' num2str(limits(4))];

X1 = X(:,1);
X2 = X(:,2);

x = X1(T);
y = X2(T);

eval(['auxmatrix = ' rule ';'])

ruledConnectivity = all(auxmatrix,2);