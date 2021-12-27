function plotHandle = plotPMLnodes(X,sigmaPos,axesHandle)

plotHandle = cell(2,1);
hold on
plotHandle{1} = plot(axesHandle,X(sigmaPos{1},1),X(sigmaPos{1},2),...
    'ro','markerFaceColor','r','markerSize',4,'Visible','off');
plotHandle{2} = plot(axesHandle,X(sigmaPos{2},1),X(sigmaPos{2},2),...
    's','markerFaceColor','none','markerSize',4,'Visible','off');
hold off