
data = guidata(handles.MainFigure);

hmax = 25;
hmin = 1;
ymax = 366;
bottom = data.bottom.value;

posmax = find(bottom > hmax);
posmin = find(bottom < hmin);

hold on
plotmax = plot(data.mesh.X(posmax,1),data.mesh.X(posmax,2),'r*');
plotmin = plot(data.mesh.X(posmin,1),data.mesh.X(posmin,2),'b*');

figure
bottom(posmax) = hmax;
bottom(posmin) = hmin;
bottom(data.mesh.X(:,2) >= ymax) = hmax;
plotSolution(data.mesh.X,data.mesh.T,bottom,data.mesh.referenceElement)

data.bottom.value = bottom;
guidata(handles.MainFigure,data)

% fid = fopen('bat_mataro[background].txt','w');
% fprintf(fid,'%i\n',numel(data.bottom.value));
% fprintf(fid,'%g %g %g\n',[data.mesh.X,data.bottom.value]');
% fclose(fid);

