function Mexpand = expandMatrix(M)

Mexpand = zeros(size(M)*2);

I = [1 0; 0 1];

for i = 1:size(M,1)
    for j = 1:size(M,2)
    Mexpand(2*i-1:2*i,2*j-1:2*j) = M(i,j)*I;
    end
end