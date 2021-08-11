% Exact.m
function [suma]=exactSolScatteringCercle(m,theta,k,a)% DGHFull(k,nDeg,Orden_NRBC)
%
% k: wave number
% m: radius of points (vector)
% theta: valores de theta para los puntos (vector)
% a: obstacle radius
% suma: analytical solution at these points

P0=1;
% switch k
%     case {1,2}
%         n_terminos=2*k+4;
%     case {3,4}
%         n_terminos=2*k+3;
%     case {5,6}
%         n_terminos=2*k+2;
%     case {7,8}
%         n_terminos=2*k+1;
%     case {9,10}
%         n_terminos=2*k-1;
%     case {11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40}
%         n_terminos=3*k/2+3;
%     otherwise 
%         if k<1
%             n_terminos=10;
%         else
%             n_terminos=60;
%         end
% end
% n_terminos=1*n_terminos;
n_terminos = 100;
nm=length(m);
suma=zeros(nm,1);
im=sqrt(-1);

for j=1:nm
    r=m(j,1);
    for nu=0:n_terminos%-n_terminos:n_terminos
        
        Ja = besselj(nu,k*a);
        Ja1 = besselj(nu + 1, a*k);
        dJa = (nu*Ja)/a - k*Ja1;
        Ya = bessely(nu,k*a);
        Ya1 = bessely(nu + 1,k*a);
        dHa = (nu*Ja)/a - k*Ya1*1i - k*Ja1 + (nu*Ya*1i)/a;

        Jr = besselj(nu,k*r);
        Yr = bessely(nu,k*r);
        Hr = Jr + im * Yr;
        
        cociente = dJa*Hr/dHa;
        
        if nu==0
            e=1;
        else
            e=2;
        end
        suma(j,1)=suma(j,1)-e*im^nu*cociente*cos(nu*theta(j,1));%exp(im*i*theta(j,1));
    end
end
suma=P0*suma;