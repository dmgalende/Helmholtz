function Sgauss = berkoffPMLSource(sigma,U)

Sgauss = zeros(size(U));

% Bonnet and Poupaud PML (PDE version)---------------------------------
Sgauss(:,1) = sigma(:,1).*U(:,1) + (sigma(:,2)-sigma(:,1)).*U(:,5);
Sgauss(:,2) = sigma(:,1).*U(:,2);
Sgauss(:,3) = sigma(:,2).*U(:,3);
Sgauss(:,4) = 0;
Sgauss(:,5) = sigma(:,2).*U(:,5);

