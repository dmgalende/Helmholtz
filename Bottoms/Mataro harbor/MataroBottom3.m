function h = MataroBottom3(X)

%MATLAB 2013 or newer required

% Real Bathymetry
xyzdata = load('MataroBottomRotated.mat');

% Bottom mesh
F = scatteredInterpolant(xyzdata.xyz(:,1),xyzdata.xyz(:,2),xyzdata.xyz(:,3),'linear','none');
h = F(X(:,1),X(:,2));

% Max bottom
% h(h > 25) = 25;

% Fix bad low values
h(h < 3) = 3;

% Find bugs (NaN values)
posnan = isnan(h);
if any(posnan)
    disp('There are NaN values to be fixed!')
end
