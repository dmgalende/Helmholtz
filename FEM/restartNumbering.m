function X_r = restartNumbering(X)

s = size(X);
n = length(X(:));
[v, ind] = sort(X(:));
values_back = 1:n;
ind_back(ind) = values_back;
newnum = zeros(n,1);
newnum(1) = 1;
p = 2;
for i = 2:n
    if v(i) == v(i-1)
        newnum(i) = newnum(i-1);
    else
        newnum(i) = p;
        p = p + 1;
    end
end
v_new = newnum(ind_back);
X_r = reshape(v_new, s);