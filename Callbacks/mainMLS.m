function approx = mainMLS(X,data,auxRho,elemSize,handle)

%_________________________________________________________________
%
% 2D functional approximation by Moving Least-Squares (MLS) 
%_________________________________________________________________

%Particle data and grid limits
zdata = data(:,3);
xydata = data(:,1:2);
xmax = max(xydata(:,1)); xmin = min(xydata(:,1)); disx = xmax - xmin;
ymax = max(xydata(:,2)); ymin = min(xydata(:,2)); disy = ymax - ymin;
coef = 0.01;
ax = xmin - coef*disx; bx = xmax + coef*disx;
ay = ymin - coef*disy; by = ymax + coef*disy;
nx = round(abs((bx-ax))/elemSize);
ny = round(abs((by-ay))/elemSize);

%Localization cells
setOutput({'Creating grid data...'},handle)
cells = particlesAtCells(xydata,nx,ny,ax,bx,ay,by);

%Dilation parameters
nOfParticles = size(data,1);
rho = auxRho*ones(nOfParticles,1); %constant rho

%MLS approximation of the data
setOutput({'Doing MLS approximation...'},handle)
nOfPoints = size(X,1);
approx = zeros(nOfPoints,1);
for i = 1:nOfPoints
    x = X(i,:);
    c = computeMLScoefficients(x,zdata,xydata,rho,cells);
    approx(i) = c'*P(x);   
end

%Check approximation
poorPos = find(isnan(approx),1);
if ~isempty(poorPos)
    setOutput({'MLS fails! Bad nodes were detected. Try increasing the rho parameter'},handle)
end











