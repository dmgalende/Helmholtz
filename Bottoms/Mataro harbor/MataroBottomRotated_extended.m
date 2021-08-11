function h = MataroBottomRotated_extended(X)

%MATLAB 2013 or newer required

% Real Bathymetry
data = load('datamataro_bottom_extendido.mat');
data = data.data;

% Bottom mesh
F = scatteredInterpolant(data.mesh.X(:,1),data.mesh.X(:,2),data.bottom.value,'linear','nearest');
h = F(X(:,1),X(:,2));

% Max bottom
h(h > 25) = 25;

% Fix bad low values
h(h < 3) = 3;

% Find bugs (NaN values)
posnan = isnan(h);
if any(posnan)
    disp('There are NaN values to be fixed!')
end
