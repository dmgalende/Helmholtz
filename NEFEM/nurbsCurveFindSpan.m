function res = nurbsCurveFindSpan(u,U)
%
% res = nurbsCurveFindSpan(u,U)
%
% Find the knot span for the parameter u
%

% NURBS info
p = length(find(U==U(1))) - 1;
m = length(U) - 1;
n = m - p - 1;

tol = 1e-5;

if abs(u-U(1)) < tol
    res = p+1;
elseif abs(u-U(end)) < tol
    res = n+1;
else    
    % Binary search
    low = p;
    high = n+1;
    mid = (low + high)/2;
    
    while (u<U(round(mid)) || u>=U(round(mid)+1))
        if u<U(round(mid))
            high = mid;
        else
            low = mid;
        end
        mid = (low + high)/2;
    end
    
    res = round(mid);
end    