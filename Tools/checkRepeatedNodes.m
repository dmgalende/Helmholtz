rep = [];
for i = 1:size(X,1)
    x = X(i,1);
    y = X(i,2);
    for j = i+1:size(X,1)
        yr = X(j,2);
        xr = X(j,1);
        if y == yr && x == xr
            rep = [rep ; i j];
        end
    end
end
    