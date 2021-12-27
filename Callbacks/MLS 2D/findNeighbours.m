function isp = findNeighbours(x,cells,Xpar,rho) 
% isp = findNeighbours(x,cells,X,rho) 
% Output: vector with the indexes of the neighbouring particles

nx = cells.Ncelx;
ny = cells.Ncely;
aux = cells.domain; 
ax = aux(1); bx=aux(2); ay=aux(3); by = aux(4);

hx = (bx-ax)/nx; hy = (by-ay)/ny;
ip = min(floor((x(1)-ax)/hx)+1,nx);
jp = min(floor((x(2)-ay)/hy)+1,ny);

rhomax = max(rho);
lx = floor(rhomax/hx)+1 ;
ly = floor(rhomax/hy)+1 ;

isp0 = [];

particel = cells.particel; nparcel=cells.nparcel;
for i=max(ip-lx,1):min(ip+lx,nx)
   for j=max(jp-ly,1):min(jp+ly,ny)
      k = (i-1)*ny+j;
      nk = nparcel(k);
      if nk>0 
         isp0 = [isp0  full(particel(k,1:nk))];  
      end
   end
end

if ~isempty(isp0)
   isp1 = isp0(abs(Xpar(isp0,1)-x(1))<rho(isp0));
else
   isp = [];
   return;
end

if ~isempty(isp1)
   isp  = isp1(abs(Xpar(isp1,2)-x(2))<rho(isp1));
else 
   isp = [];
end