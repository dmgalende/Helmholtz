function FluxU_PEC = berkoffFluxU_PEC(UgL,x,n,t,c,alpha_face,incidentWave,ABS)

if size(n,1)==1
    n = repmat(n, size(UgL,1), 1);
end
n1 = n(:,1);
n2 = n(:,2);

sqrt_alpha_face = sqrt(alpha_face);
FluxU_PEC = zeros(size(UgL));
[phix,phiy,phit] = berkoffIncidentWave(t,x,incidentWave);
g = -(n1.*phix + n2.*phiy +ABS./c.*phit);
Delta = (n1.*UgL(:,2) + n2.*UgL(:,3));
U1_mas = (UgL(:,1)+sqrt_alpha_face.*(g-Delta))./(1+sqrt_alpha_face*ABS./c);
Delta_mas = (g + ABS./c.*(sqrt_alpha_face.*Delta-UgL(:,1)))./(1+ABS./c.*sqrt_alpha_face);

% Left
FluxU_PEC(:,1) = 0.5* sqrt_alpha_face.*((UgL(:,1) - U1_mas) - sqrt_alpha_face.*(Delta + Delta_mas));
FluxU_PEC(:,2) = 0.5*n1.*(-(UgL(:,1) + U1_mas) + sqrt_alpha_face.*(Delta - Delta_mas));
FluxU_PEC(:,3) = 0.5*n2.*(-(UgL(:,1) + U1_mas) + sqrt_alpha_face.*(Delta - Delta_mas));
