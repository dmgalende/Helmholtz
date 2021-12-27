function res = phi1dfun(z) 
 r = abs(z); 
 if (r>=1)  
    res = 0; 
 elseif (r>=0.5) 
    res = ((-r+1)^3)*4/3;  
 else 
    res = 2/3-4*r^2*(1-r); 
 end 
 
%Funcion Gaussiana 
%res = (exp(-9*((r/rho)^2))-exp(-9))/(1-exp(-9)); 
