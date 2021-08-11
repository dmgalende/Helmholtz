function [phix,phiy,phit] = berkoffIncidentWave(t,coordFace,incidentWave)
%
%
% phi = A cos(omega*t - kVectorX)
%
%

k = incidentWave.waveNumberValue;
A = incidentWave.amplitude;
omega = 2*pi/incidentWave.period;
kx = k * cos(pi/180*incidentWave.direction);
ky = k * sin(pi/180*incidentWave.direction);
kVectorX = coordFace(:,1)*kx + coordFace(:,2)*ky;

nOfPoints = length(coordFace(:,1));
phix = zeros(nOfPoints,1);
phiy = zeros(nOfPoints,1);
phit = zeros(nOfPoints,1);

% modulation
aux = omega*t/(2*pi);
if aux<=1
    modulation = aux;
else
    modulation = 1;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iPoint = 1:nOfPoints
    aux = omega*t - kVectorX(iPoint);
%     if aux <= 0
%         phix(iPoint) = 0;
%         phiy(iPoint) = 0;
%     elseif aux <= 2*pi
%         phix(iPoint) = kx * A * (aux/(2*pi))*sin(aux);
%         phiy(iPoint) = ky * A *(aux/(2*pi))*sin(aux);
%     else
%         phix(iPoint) = kx * A *sin(aux);
%         phiy(iPoint) = ky * A *sin(aux);
%     end

    phix(iPoint) = modulation * kx * A *sin(aux);
    phiy(iPoint) = modulation * ky * A *sin(aux);
    phit(iPoint) = modulation * (-omega * A *sin(aux));

end
