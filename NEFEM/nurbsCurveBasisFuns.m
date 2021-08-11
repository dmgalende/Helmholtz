function N = nurbsCurveBasisFuns(i,u,U)
%
% N = nurbsCurveBasisFuns(i,u,U)
% 
% Non-zero basis functions for a NURBS curve
%

p = length(find(U==U(1))) - 1;

N = zeros(1,p+1);
left = zeros(1,p+1);
right = zeros(1,p+1);

N(1)=1;

for j=1:p
    left(j) = u - U(i+1-j);
    right(j) = U(i+j) - u;
    res = 0;
    
    for k=1:j
        aux = N(k)/(right(k) + left(j-k+1));
        N(k) = res + right(k)*aux;
        res = left(j-k+1)*aux;
    end
    
    N(j+1) = res;
end
