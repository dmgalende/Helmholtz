function h = MataroBottom2(X)

% Real Bathymetry
xyzdata = load('MataroBottomRotated.mat');

% Bottom mesh
h = griddata(xyzdata.xyz(:,1),xyzdata.xyz(:,2),xyzdata.xyz(:,3),X(:,1),X(:,2),'cubic');

% Max bottom
h(h > 25) = 25;

% Fix bad low values
h(h < 3) = 3;

% Find bugs (NaN values)
posnan = isnan(h);
if any(posnan)
    disp('There are NaN values to be fixed!')
end
