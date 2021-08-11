function cells = particlesAtCells(Xpar,nx,ny,ax,bx,ay,by)

 n = size(Xpar,1);
 
 particel = spalloc((nx+1)*(ny+1),n,n);
 nparcel  = zeros((nx+1)*(ny+1),1);
 
 for k=1:n
     ncel = nOfCell(Xpar(k,:),nx,ny,ax,bx,ay,by);
     nparcel(ncel) = nparcel(ncel)+1;
     particel(ncel,nparcel(ncel)) = k;
 end
 
 cells.nparcel = nparcel;
 cells.particel = particel;
 cells.Ncelx = nx;
 cells.Ncely = ny;
 cells.domain = [ax,bx,ay,by];
 
%______________________________________________________
function nc = nOfCell(x,nx,ny,ax,bx,ay,by)

 hx = (bx-ax)/nx; hy = (by-ay)/ny;
 i = min(floor((x(1)-ax)/hx)+1,nx);
 j = min(floor((x(2)-ay)/hy)+1,ny);
 nc = (i-1)*ny+j;
 