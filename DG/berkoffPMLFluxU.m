function [FluxUL,FluxUR] = berkoffPMLFluxU(uL, uR, n, alpha_face)

if size(n,1)==1
    n = repmat(n, size(uL,1), 1);
end
n1 = n(:,1);
n2 = n(:,2);

% Initialization
FluxUL = zeros(size(uL));
FluxUR = FluxUL;

Delta = n1.*uL(:,2) +n2.*uL(:,3);
Delta_mas = n1.*uR(:,2) +n2.*uR(:,3);
U3mean = 0.5*alpha_face.*( uL(:,3) + uR(:,3) );
sqrt_alpha_face = sqrt(alpha_face);

% Left
FluxUL(:,1) = 0.5* sqrt_alpha_face.*((uL(:,1) - uR(:,1)) - sqrt_alpha_face.*(Delta + Delta_mas));
FluxUL(:,2) = 0.5*n1.*(-(uL(:,1) + uR(:,1)) + sqrt_alpha_face.*(Delta - Delta_mas));
FluxUL(:,3) = 0.5*n2.*(-(uL(:,1) + uR(:,1)) + sqrt_alpha_face.*(Delta - Delta_mas));
FluxUL(:,5) = -n2.*U3mean;

n1 = -n1;
n2 = -n2;
Delta = n1.*uR(:,2) +n2.*uR(:,3);
Delta_mas = n1.*uL(:,2) +n2.*uL(:,3);

% Right
FluxUR(:,1) = 0.5* sqrt_alpha_face.*((uR(:,1) - uL(:,1)) - sqrt_alpha_face.*(Delta + Delta_mas));
FluxUR(:,2) = 0.5*n1.*(-(uR(:,1) + uL(:,1)) + sqrt_alpha_face.*(Delta - Delta_mas));
FluxUR(:,3) = 0.5*n2.*(-(uR(:,1) + uL(:,1)) + sqrt_alpha_face.*(Delta - Delta_mas));
FluxUR(:,5) = -n2.*U3mean;



