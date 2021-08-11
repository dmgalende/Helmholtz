function plotPmap(mesh,p_adaptivity)

% Mesh info
nameNodal = mesh.fieldNames{mesh.indexElemPosCon(2)};
nameCon = mesh.fieldNames{mesh.indexElemPosCon(3)};
X = mesh.(nameNodal);
T = mesh.(nameCon);
nOfElements = size(T,1);
elem_p = p_adaptivity.elem_p;
different_p = p_adaptivity.different_p;
leg = cell(size(different_p));
hand = zeros(size(different_p));

% create legend
for i = 1:numel(leg)
    leg{i} = num2str(different_p(i));
end

% Loop in elements
for ielem = 1:nOfElements

    % Interpolate solution and position at interpolation points
    Te = T(ielem,1:3);
    Xe = X(Te,:);
    u = elem_p(ielem);
    % Plot interpolated solution in the element
    hold on
    hand(different_p == u) = patch(Xe(:,1),Xe(:,2),u);
    hold off
end
axis equal
legend(hand,leg)
end
