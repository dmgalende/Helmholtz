function h = bottomFunctionShoaling_provaPuerto(X)

x = X(:,1);
y = X(:,2);
longitud = max(x)-min(x);
pml = 113;
h = 5 + (x-min(x))/longitud*25;
h(X(:,1)<=min(x)+pml | X(:,1)>=max(x)-pml | X(:,2)>=max(y)-pml) = 30;


