function h = BarcelonaBottom(X)

bottomData = './../../Saves/Barcelona harbor/BarT6SE_h15-30_PMLn2R20_DG_P1.mat';

load(bottomData)
D = scatteredInterpolant(data.mesh.X(:,1),data.mesh.X(:,2),data.bottom.value,'linear','nearest');
h = D(X(:,1),X(:,2));
if any(h<0), warning('Negative values of bottom interpolation detected'), h(h<0) = min(h(h>0)); end