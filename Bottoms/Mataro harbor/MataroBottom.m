function h = MataroBottom(X)

% Real Bathymetry
xyzdata = load('MataroBottom.mat');

% Bottom mesh
h = griddata(xyzdata.xyz(:,1),xyzdata.xyz(:,2),xyzdata.xyz(:,3),X(:,1),X(:,2),'cubic');

% Max bottom
h(h > 25) = 25;

% Fix bad low values
h(h < 3) = 3;

% Find bugs (NaN values)
posnan = isnan(h);
if any(posnan)
    lims = [-640 -625   160 170];
    pos = find(posnan & X(:,2) > lims(3) & X(:,2) < lims(4) & X(:,1) > lims(1) & X(:,1) < lims(2));
    h(pos) = 3;
    posnan(pos) = false;
    lims = [ -450 -415 140 146];
    pos = find(posnan & X(:,2) > lims(3) & X(:,2) < lims(4) & X(:,1) > lims(1) & X(:,1) < lims(2));
    h(pos) = 3;
    posnan(pos) = false;
    lims = [ 295 300 -329 -328];
    pos = find(posnan & X(:,2) > lims(3) & X(:,2) < lims(4) & X(:,1) > lims(1) & X(:,1) < lims(2));
    h(pos) = 3;
    posnan(pos) = false;
    if any(posnan)
        disp('There are NaN values to be fixed!')
    end
end
