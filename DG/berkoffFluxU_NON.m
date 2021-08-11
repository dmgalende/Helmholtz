function FluxU_NON = berkoffFluxU_NON(UgL,n,alpha_face)

if size(n,1)==1
    n = repmat(n, size(UgL,1), 1);
end
n1 = n(:,1);
n2 = n(:,2);

sqrt_alpha_face = sqrt(alpha_face);
FluxU_NON = zeros(size(UgL));

Delta = (n1.*UgL(:,2) + n2.*UgL(:,3));
U1_mas = UgL(:,1) - sqrt_alpha_face .* Delta;
Delta_mas = 0;

% Left
FluxU_NON(:,1) = 0.5* sqrt_alpha_face.*((UgL(:,1) - U1_mas) - sqrt_alpha_face.*(Delta + Delta_mas));
FluxU_NON(:,2) = 0.5*n1.*(-(UgL(:,1) + U1_mas) + sqrt_alpha_face.*(Delta - Delta_mas));
FluxU_NON(:,3) = 0.5*n2.*(-(UgL(:,1) + U1_mas) + sqrt_alpha_face.*(Delta - Delta_mas));
