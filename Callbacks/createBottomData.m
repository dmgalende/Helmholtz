function  flag = createBottomData(file,ext,handles)

handle = handles.run_wipOutput;
data = guidata(handle);

if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
    nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
elseif strcmp(data.computation,'NEFEM')
    nameNodal = 'X';
end

switch ext
    case '.mat'
        auxBottom = load(file);
        auxBottomField = fieldnames(auxBottom);
        bottom = auxBottom.(auxBottomField{:});
        if numel(auxBottom) ~= 1
            msje = {['??? The .mat file has to contain only the bottom variable. '...
                'Check the Berkhoff GUI help for details']};
            setOutput(msje,handle), error('Error in Berkhoff GUI')
        elseif ~isnumeric(bottom) || all(size(bottom) ~= ...
                size(data.mesh.(nameNodal),1)) || all(size(bottom) ~= 1)
            msje ={['??? The bottom cant be read. Make sure it coincides with the '...
                'selected mesh']};
            setOutput(msje,handle), error('Error in Berkhoff GUI')
        end
    case '.m'
        bottomFun = str2func(file(1:end-2));
        if any([nargin(bottomFun) nargout(bottomFun)] ~= 1)
            msje = {['??? The selected m-file has to have only one input/output '...
                'argument. Check the Berkhoff GUI help for details']};
            setOutput(msje,handle), error('Error in Berkhoff GUI')
        end
        bottom = bottomFun(data.mesh.(nameNodal));
        if ~isnumeric(bottom) || all(size(bottom) ~= ...
                size(data.mesh.(nameNodal),1)) || all(size(bottom) ~= 1)
            msje = {['??? The computed bottom doesnt coincide with the selected '...
                'mesh. Check the Berkhoff GUI help for details']};
            setOutput(msje,handle), error('Error in Berkhoff GUI')
        end
    case '.txt'
        setOutput({'Reading file...'},handle)
        try
            fid = fopen(file,'r');
            nOfPoints = str2num(fgetl(fid));
            if isscalar(nOfPoints)
                auxBottom = zeros(nOfPoints,3);
                for ipoint = 1:nOfPoints
                    auxBottom(ipoint,:) = str2num(fgetl(fid));
                end
            elseif numel(nOfPoints) == 3
                auxBottom(1,:) = nOfPoints;
                ipoint = 2;
                line = fgetl(fid);
                while line ~= -1
                    auxBottom(ipoint,:) = str2num(line);
                    line = fgetl(fid);
                    ipoint = ipoint + 1;
                end
            else
                msje = {['??? Bottom data cant be read. The first row should '...
                         'be the total number of points or the (x,y,depth) data '...
                         'of the first one']};
                setOutput(msje,handle), error('Error in Berkhoff GUI')
            end
            fclose(fid);
        catch
            msje = {['??? Bottom data cant be read. Make sure that each row '...
                     'of the selected file contains the (x,y,depth) data of '...
                     'the point']};
            setOutput(msje,handle), error('Error in Berkhoff GUI')
        end
        handleData = plotData(auxBottom,handles.axesHandle);
        figBottom = plot3Bottom(auxBottom);
        setOutput({'Enter MLS parameters'},handle)
        [rho,elemSize] = inputfun('Rho parameter: ','Grid size: ','MLS options',handle);
        delete(handleData)
        try close(figBottom), end
        if isempty([rho elemSize])
            flag = 0;
            return
        end
        bottom = mainMLS(data.mesh.(nameNodal),auxBottom,rho,elemSize,handle);
    otherwise
        msje = {['??? File extension unrecognized. Check the Berkhoff '...
                 'GUI help for details']};
        setOutput(msje,handle), error('Error in Berkhoff GUI')
end

if isfield(data,'bottom')
    data = rmfield(data,'bottom');
end
data.bottom.value = bottom;

data.constantBottomFlag = false;

guidata(handle,data);

flag = 1;


%_____________________________________
function [rho,elemSize] = inputfun(title1,title2,title,handle)

default = {'none' 'none'};
options.WindowStyle = 'normal';
answer = inputdlg({title1 title2},title,1,default,options);
if isempty(answer)
    rho = [];
    elemSize = [];
    return
end

rho = str2num(answer{1});
elemSize = str2num(answer{2});
if any([isempty(rho) isempty(elemSize) ~isscalar(rho) ~isscalar(elemSize)])
    setOutput({'??? Incorrect parameters'},handle)
    [rho,elemSize] = inputfun(title1,title2,title,handle);
end


%_____________________________________
function handleData = plotData(data,handle)

hold on
handleData = plot(handle,data(:,1),data(:,2),'o','markerFaceColor','r','markerSize',3);
hold off


%_____________________________________
function figBottom = plot3Bottom(data)

figBottom = figure('Name','Bottom'); 
plot3(data(:,1),data(:,2),-data(:,3),'o');
xlabel('x'), ylabel('y'), zlabel('Depth')



