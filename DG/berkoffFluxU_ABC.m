function FluxU_ABC = berkoffFluxU_ABC(UgL,n,c,alpha_face,R)

if size(n,1)==1
    n = repmat(n, size(UgL,1), 1);
end
n1 = n(:,1);
n2 = n(:,2);

sqrt_alpha_face = sqrt(alpha_face);
FluxU_ABC = zeros(size(UgL));

g = 0;
Delta = (n1.*UgL(:,2) + n2.*UgL(:,3));
U1_mas = (UgL(:,1)+sqrt_alpha_face.*(g-Delta-UgL(:,4)/2/R))./(1+sqrt_alpha_face*1./c);
Delta_mas = (g + 1./c.*(sqrt_alpha_face.*Delta+sqrt_alpha_face.*UgL(:,4)/2/R-UgL(:,1)))./(1+1./c.*sqrt_alpha_face);

% Left
FluxU_ABC(:,1) = 0.5* sqrt_alpha_face.*((UgL(:,1) - U1_mas) - sqrt_alpha_face.*(Delta + Delta_mas));
FluxU_ABC(:,2) = 0.5*n1.*(-(UgL(:,1) + U1_mas) + sqrt_alpha_face.*(Delta - Delta_mas));
FluxU_ABC(:,3) = 0.5*n2.*(-(UgL(:,1) + U1_mas) + sqrt_alpha_face.*(Delta - Delta_mas));
