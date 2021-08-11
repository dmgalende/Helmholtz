function [Cxi,Ceta,Czeta] = ...
    computeReferenceElementConvectionMatrices(nDeg,invV)
%
% [Cxi,Ceta,Czeta] = computeReferenceElementConvectionMatrices(nDeg,invV)
%
% Function for the computation of the convection matrices in
% the reference element, with the ortogonal polinomial base
% of degree less or equal to nDeg.
%
% Input:
% nDeg:  degree of interpolation
% invV:  inverse of the vandermonde matrix
%
% Output:
% [Cxi,Ceta,Czeta]: convection matrices in local coordinates
%                   size is nOfNodes X nOfNodes
%

switch nargout
    case 1,
        Cxi = compute_Cs1D(nDeg,invV);
    case 2,
        [Cxi,Ceta] = compute_Cs2D(nDeg,invV);
    case 3,
        [Cxi,Ceta,Czeta] = compute_Cs3D(nDeg,invV);
    otherwise
        error('Wrong number of output arguments in computeReferenceElementConvectionMatrices')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cxi = compute_Cs1D(nDeg,invV)

nOfNodes = nDeg+1;%number of nodes/polynomials
Cxi  = zeros(nOfNodes);

[z,w] = gaussLegendre(nDeg+1,-1,1);

nIP = length(w); %number of integration points
%Integration over [-1,1]
for i = 1:nIP
    x = z(i);
    weight = w(i);
    [p,p_xi] = orthopoly1D_deriv(x,nDeg);
    Cxi  = Cxi  + weight*p*p_xi';
end

Cxi  = (invV')*Cxi*invV;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cxi,Ceta] = compute_Cs2D(nDeg,invV)

nOfNodes = (nDeg+1)*(nDeg+2)/2 ;%number of nodes/polynomials
Cxi  = zeros(nOfNodes);
Ceta = zeros(nOfNodes);

[z,w] = gaussLegendre(nDeg+1,-1,1);

nIP = length(w); %number of integration points in each direction
%Integration over [-1,1]^2
for i = 1:nIP
    for j = 1:nIP
        x = [z(i),z(j)]; %(r,s)
        weight = (w(i)*w(j))*(1-x(2))/2;
        [p,p_xi,p_eta] = orthopoly2D_deriv_rst(x,nDeg);
        Cxi  = Cxi  + weight*p*p_xi';
        Ceta = Ceta + weight*p*p_eta';
    end
end

Cxi  = (invV')*Cxi*invV;
Ceta = (invV')*Ceta*invV;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cxi,Ceta,Czet] = compute_Cs3D(nDeg,invV)

nOfNodes = (nDeg+1)*(nDeg+2)*(nDeg+3)/6 ;%number of nodes/polynomials
Cxi  = zeros(nOfNodes);
Ceta = zeros(nOfNodes);
Czet = zeros(nOfNodes);

[z,w] = gaussLegendre(nDeg+1,-1,1);

nIP = length(w); %number of integration points in each direction
%Integration over [-1,1]^3
for i = 1:nIP
    for j = 1:nIP
        for k = 1:nIP
            x = [z(i),z(j),z(k)]; %(r,s,t)
            weight = (w(i)*w(j)*w(k))*jacobian(x);
            [p,p_xi,p_eta,p_zet] = orthopoly3D_deriv_rst(x,nDeg);
            Cxi  = Cxi  + weight*p*p_xi';
            Ceta = Ceta + weight*p*p_eta';
            Czet = Czet + weight*p*p_zet';
        end
    end
end

Cxi  = (invV')*Cxi*invV;
Ceta = (invV')*Ceta*invV;
Czet = (invV')*Czet*invV;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = jacobian(x)
r = x(1); s = x(2); t = x(3);
J = ((1-s)/2)*((1-t)/2)^2;
