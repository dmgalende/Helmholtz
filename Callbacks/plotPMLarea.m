function plotHandle = plotPMLarea(boundaries,lon,axesHandle)

%Empty boundaries
pos = [];
sum = 0;
for i = 1:4
    if ~isempty(boundaries{i})
        pos = [pos i];
        sum = sum + size(boundaries{i},1);
    end
end

%Area plot
PMLlon = [lon,-lon,-lon,lon];
coord = [1 2 1 2];
plotHandle = zeros(sum,1);
u = 1;
for i = pos
    for k = 1:size(boundaries{i},1)
        x = boundaries{i}(k,1:2)';
        y = boundaries{i}(k,3);
        dy = y + PMLlon(i);
        area2plot = [[x; x(2); x(1)] [y; y; dy; dy]];
        tcolor(1,1,1:3) = [1 0 0];
        plotHandle(u) = patch(area2plot(:,coord(i)),area2plot(:,3-coord(i)),...
            tcolor,'parent',axesHandle,'EdgeColor','none');
        u = u + 1;
    end
end