function h = BarcelonaBottomT6_trasl(X)

% Real Bathymetry
xyzdata = load('BarcelonaBat.mat');
xyzdata.xyz(:,1) = xyzdata.xyz(:,1) - 424096;
xyzdata.xyz(:,2) = xyzdata.xyz(:,2) - 4568110;
% Rotate
RotPoint = [429591.7298481572 - 424096,4569045.335042576 - 4568110];
RotAngle = (180 - 67.4056)*pi/180;
R = [cos(RotAngle) -sin(RotAngle)
     sin(RotAngle) cos(RotAngle)];
oneVec = ones(1,size(xyzdata.xyz,1));
c = [RotPoint(1)*oneVec ; RotPoint(2)*oneVec];
xy = (R*(xyzdata.xyz(:,1:2)' - c)) + c;

% Bottom mesh
h = griddata(xy(1,:),xy(2,:),-xyzdata.xyz(:,3),X(:,1),X(:,2),'cubic');

% Max bottom
h(h > 30) = 30;

% Fix bad low values
h(h < 3) = 3;

% Find bugs (NaN values)
posnan = isnan(h);
if any(posnan)
    lims = [ 0.4172-0.424096    0.4173-0.424096    4.5663-4.568110    4.5664-4.568110]*1e6;
    pos = find(posnan & X(:,2) > lims(3) & X(:,2) < lims(4) & X(:,1) > lims(1) & X(:,1) < lims(2));
    h(pos) = 10;
    posnan(pos) = false;
    lims = [0.4173-0.424096    0.4174-0.424096    4.5662-4.568110    4.5663-4.568110]*1e6;
    pos = find(posnan & X(:,2) > lims(3) & X(:,2) < lims(4) & X(:,1) > lims(1) & X(:,1) < lims(2));
    h(pos) = 7;
    posnan(pos) = false;
    if any(posnan)
        disp('There are NaN values to be fixed!')
    end
end