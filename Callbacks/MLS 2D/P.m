function res = P(x0) 
 n = size(x0,1);
 x = x0(:,1)';
 y = x0(:,2)';
 res = [ones(1,n) 
        x 
        y 
%         x.^2
%         x.*y
%         y.^2
]; 
