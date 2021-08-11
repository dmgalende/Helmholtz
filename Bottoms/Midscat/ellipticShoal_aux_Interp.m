function h = ellipticShoal_aux_Interp(X)

data = load('bottom_ref_2.mat');

D = scatteredInterpolant(data.X(:,1),data.X(:,2),data.bottom,'linear','none');
h = D(X(:,1),X(:,2));
