function c = computeMLScoefficients(x,u,X,rho,cells)

%Neigbouring particles
isp = findNeighbours(x,cells,X,rho); %indexes
Xr = X(isp,:); %coordinates
rhor = rho(isp);

%Gramm matrix for the least-squares local fitting
M = GrammMatrix(x,Xr,rho);

%Right-hand side: <P,u>
ur = u(isp); %value of u at neighbouring particles
Pu = zeros(size(M,1),1); %memory allocation
for i=1:size(Xr,1)
  Pu = Pu + P(Xr(i,:))*ur(i)*phifun(Xr(i,:)-x,rhor(i));   
end

c = M\Pu;

%______________________________________________
function M = GrammMatrix(x,Xr,rhor)
% M = GrammMatrix(x,Xr,rhor)
% Input:
%    x: point for the evaluation
%    Xr: coordinates of the neigbouring particles
%    rhor: dilation parameter of these particles
% Output:
%    M: Gramm matrix for the local Least-Squares fitting

 n = length(P([0,0])); %dimension of the polynomial space
 M = zeros(n,n); %memory allocation (full matrix)
 
 %loop for neighbouring particles
 for j=1:size(Xr,1)
   xj = Xr(j,:); %particle coordinate
   Pj = P(xj); %polinomial base at particle j
   M = M + Pj*((Pj')*phifun(x-xj,rhor(j))); %contribution to matrix
 end

 