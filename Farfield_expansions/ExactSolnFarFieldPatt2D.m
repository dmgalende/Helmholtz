
% SCATTERED FIELD AND FAR FIELD PATTERN FOR THE SCATTERED PRESSURE FROM A
% SOFT INFINITE CIRCULAR CYLINDER


function [uscex,A0ex] = ExactSolnFarFieldPatt2D(NFS, r, theta, k, r0)
% r:   Radial variable
% theta:   Angular variable
% k:    Wavenumber
% ro:   Radius of circle (Inner boundary)
% NFS:  Number of terms used in the analytical formula for the Fourier
%       Series.

% global r;
% global theta;
% global k;
% global r0;
% global N;
% global m;

N = length(r);
m = length(theta);

kro=k*r0;
imag = sqrt(-1);
Hn_ro = zeros(NFS);
Jn_ro = zeros(NFS);
cos_theta = zeros(NFS,m);
Hn_r = zeros(NFS,N);

%fprintf('\nComputing Exact Solution Metrics\n');
 for n = 0:NFS
    Hn_ro(n+1) = besselh(n,kro);
    Jn_ro(n+1) = besselj(n,kro);
    for j = 1:m
        cos_theta(n+1,j) = cos(n*theta(j));
    end
    for i=1:N
        Hn_r(n+1,i) = besselh(n,k*r(i));
    end
 end
 
% Computing the scattered field
uscex = zeros(N,m);
for i=1:N     
    for j=1:m
        sum = -(Jn_ro(1)/Hn_ro(1))*Hn_r(1,i);
        for n = 1:NFS
            sum = sum - 2*imag^n*(Jn_ro(n+1)/Hn_ro(n+1))*...
         Hn_r(n+1,i)*cos_theta(n+1,j);
        end
        uscex(i,j) = sum;
    end
end

% Computing the Far Field Pattern from Exact Solution
A0ex = zeros(m,1); 
 for j =1:m
     sum = -(Jn_ro(1)/Hn_ro(1));
     for n=1:NFS
        sum = sum - 2*Jn_ro(n+1)/Hn_ro(n+1)*cos_theta((n+1),j);
     end
     A0ex(j) = abs(sum);
 end
 A0ex = (2/(k*pi))^(1/2)*exp(-imag*pi/4)*A0ex;

%Wrap around
uscex(:,m+1) = uscex(:,1);
A0ex(m+1) = A0ex(1);
end
