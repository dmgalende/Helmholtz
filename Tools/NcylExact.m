 function [u_B,far_ex] = NcylExact(k,z,phi,r,center,M,Nfar,Bx,By)
%--------------------------------------------------------------------------
% About
%--------------------------------------------------------------------------
%
% By: Jacob Badger and Vianey Villamizar (7/4/2018)
%
% ----------------------------- Notes -------------------------------------
% This calculates the exact solution for multiple scattering by N
% cylinders. Solution taken from "Multiple Scattering" by Paul Martin
% Specifically, this function calculates the exact solution on the
% artifical boundary of obstacles, and the exact FarField expansion.
% 
% ----------------------------- Inputs ------------------------------------
% k      (Int): Wave Number
% z      (Float in interval [0,1]): Acoustic hardness of obstacles
% phi    (Int): Angle of incidence
% r      (Float): Radius of each cylinder
% center ([n_obstacle x 2] Array): (x,y) location of center of each cylinder
% M      (Int): Number of terms in calculating exact solution
% Nfar   (Int): Number of equally spaced points on which FFP is evaluated 
% Bx     ({n_obstacle x 1} Cell Array): x coordinate of each Art Bndry
%        ordered along the theta dimension (ascend order). 
% By     ({n_obstacle x 1} Cell Array): y coordinate of each Art Bndry
%        ordered along the theta dimension (ascend order).
%
% ----------------------------- Outputs -----------------------------------
% u_B    ({n_obstacle x 1} Cell Array): Value of scattered field at Art Bndry
%                                       of each cylinder
% far_ex ([Nfar x 1] Array): Farfield pattern at Nfar equally spaced pts


% ---------------------------------------------------------------------
% Input Data
% ---------------------------------------------------------------------
Xc = center(:,1);
Yc = center(:,2);
N = length(Xc);

% Calculate position of cylinders relative to one another
b = zeros(N);
beta = zeros(N);
for q = 1:N
    for el = 1:N
        if el ~= q
            xcord=Xc(el) - Xc(q);
            ycord=Yc(el) - Yc(q);
            b(el,q) = sqrt(xcord^2 + ycord^2);    % Distance between centers
            beta(el,q) = atan2(ycord, xcord);         % Relative polar angle
        end
    end
end
% beta=beta + 2*pi*(beta<0);                % Shift to interval (0,2*pi]
A = zeros(N*(2*M + 1));                          % Main Coeff Matrix
RHS = zeros(N*(2*M + 1),1);                      % Right-hand-side vector 

% ---------------------------------------------------------------------
% Construct Coeff Matrix and RHS vector
% ---------------------------------------------------------------------
for q = 1:N
    m = -M;
    for j = 1:(2*M+1)
        row1 = (q-1)*(2*M + 1) + j;
        col1 = row1;
        A(row1,col1) = z*k*(m/(k*r(q))*besselh(m,k*r(q)) - besselh(m+1,k*r(q))) + (1-z)*besselh(m,k*r(q));
        Besselj_m = besselj(m,k*r(q));
        Besselj_m1 = besselj(m+1,k*r(q));
        temp = (z*k*(m/(k*r(q))*Besselj_m - Besselj_m1) + (1-z)*Besselj_m);
        
%       If object not at origin, shift plane wave source
        if ( Xc(q) ~= 0 || Yc(q) ~= 0)
            n = -M;
            for p = 1:(2*M + 1)
                RHS(row1) = RHS(row1) ...
                    + (-(1i)^n*exp(-1i*n*phi))*temp...
                    * besselj((n-m),k*sqrt(Xc(q)^2 + Yc(q)^2))*exp(1i*(n-m)*(atan2(Yc(q),Xc(q))));
                n = n + 1;
            end
        else
%       Otherwise store source term
            RHS(row1) = (-1i^m*exp(-1i*m*phi))*temp;
        end
        
        for el = 1:N
            if el ~= q
                n = -M;
                for p = 1:(2*M + 1)
                    col1 = (el-1)*(2*M + 1) + p;
                    A(row1,col1) = (z*k*(m/(k*r(q))*Besselj_m - Besselj_m1) + (1-z)*Besselj_m)*...
                        besselh(n-m,k*b(q,el))*exp(1i*(n-m)*beta(q,el));
                    n = n + 1;
                end
            end
        end
        m = m + 1;
    end    
end

% ---------------------------------------------------------------------
% Solving System
% ---------------------------------------------------------------------

temp = A\RHS;

c = zeros(2*M + 1,N);
for q = 1:N
    c(:,q) = temp(((q-1)*(2*M+1)+1):((q)*(2*M+1)));
end

% ---------------------------------------------------------------------
% Computing Exact Solution on Artificial Boundary
% ---------------------------------------------------------------------
% Define polar coordinates w.r.t. each circle
[Nb,~] = size(Bx);       % Npts is number of points to approx on each boundary
R = cell(N,Nb);          % Nb is number of boundaries given (doesn't have to be the same as N)
theta = cell(N,Nb);
for q = 1:N
    for j = 1:Nb
        [Npts,~] = size(Bx{j});
        R{q,j} = zeros(Npts,1);
        theta{q,j} = zeros(Npts,1);
        
        % Polar coordinates of j'th boundary pts w.r.t the q'th obstacle
        R{q,j}(:) = sqrt((Bx{j}(:) - Xc(q)).^2 + (By{j}(:) - Yc(q)).^2);  
        theta{q,j}(:) = atan2(By{j}(:) - Yc(q), Bx{j}(:) - Xc(q));       
%         theta{q,j} = theta{q,j} + 2*pi*(theta{q,j}<=0);
    end
end

% Calculate exact artificial boundary
u_B = cell(Nb,1);
for q = 1:Nb                % Calculate points on boundary of qth obstacle
    [Npts,~] = size(Bx{q});
    u_B{q} = zeros(Npts,1);
    m = -M; 
    for j = 1:2*M+1         % Sum over fourier coefficients
        for p = 1:N         % Sum over obstacles
            u_B{q}(:) = u_B{q}(:) + c(j,p)*besselh(m,k*R{p,q}(:)).*exp(1i*m*theta{p,q}(:));
        end
        m = m + 1;
    end
end

% ---------------------------------------------------------------------
% Computing Far-field pattern 
% ---------------------------------------------------------------------
theta_far = linspace(2*pi,0,Nfar).';                % Grid to calc farfield on
s=0;
m=-M;
for p=1:2*M+1
    for q = 1:N
        s = s + (-1i)^m*exp(1i*m*theta_far)...
            .*(c(p,q).*exp(-1i*k*(Xc(q)*cos(theta_far) + Yc(q)*sin(theta_far))));
    end
    m = m+1;
end
far_ex = (1-1i)/sqrt(pi) * s;
