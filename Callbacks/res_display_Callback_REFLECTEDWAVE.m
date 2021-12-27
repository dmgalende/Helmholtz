function res_display_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

tags = {
        'realPhir_handle'
        'imPhir_handle'
        'absPhir_handle'
        'realPhi_handle'
        'imPhi_handle'
        'absPhi_handle'
        'HH0_handle'
        'realPhi0'
        };

% Results figure
displayResults = figure('Name','Results','MenuBar','none','Unit',...
    'normalized','Position',[.35 .4 .25 .35],'WindowStyle','modal');

% Visualitation method 
visualitation_panel = uipanel(displayResults,'Title','Visualitation method',...
    'Unit','normalized','Position',[0 .8 1 .19]);
interpCheck_handle = uicontrol(visualitation_panel,'Style','checkbox',...
    'Unit','normalized','String','Interpolate solution','Position',[.05 .35 .47 .3],...
    'Callback',@interpCheck_Callback);
interpDeg_handle = uicontrol(visualitation_panel,'Style','edit',...
    'Unit','normalized','String','20','Position',[.55 .35 .08 .3]);
interpDegText_handle = uicontrol(visualitation_panel,'Style','text',...
    'Unit','normalized','String','Element degree','Position',[.64 .32 .32 .3]);
if any(strcmp(data.computation,{'CDG' 'DG'}))
    interpCheck = true;
    set(interpCheck_handle,'Value',1,'Enable','off')
else
    interpCheck = false;
    set(interpCheck_handle,'Value',0)
    set([interpDeg_handle interpDegText_handle],'Enable','off')
end

% Options for reflected potential
phir_panel = uipanel(displayResults,'Title','Reflected potential',...
    'Unit','normalized','Position',[0 .6 1 .19]);
realPhir_handle = uicontrol(phir_panel,'Style','checkbox','Unit','normalized',...
    'String','real','Position',[.05 .35 .2 .3],'Value',1);
imPhir_handle = uicontrol(phir_panel,'Style','checkbox','Unit','normalized',...
    'String','imag','Position',[.4 .35 .2 .3]);
absPhir_handle = uicontrol(phir_panel,'Style','checkbox','Unit','normalized',...
    'String','abs','Position',[.78 .35 .2 .3]);

% Options for total potential
phi_panel = uipanel(displayResults,'Title','Total potential',...
    'Unit','normalized','Position',[0 .4 1 .19]);
realPhi_handle = uicontrol(phi_panel,'Style','checkbox','Unit','normalized',...
    'String','real','Position',[.05 .35 .2 .3],'Value',1);
imPhi_handle = uicontrol(phi_panel,'Style','checkbox','Unit','normalized',...
    'String','imag','Position',[.4 .35 .2 .3]);
absPhi_handle = uicontrol(phi_panel,'Style','checkbox','Unit','normalized',...
    'String','abs','Position',[.78 .35 .2 .3]);

% Various options 
various_panel = uipanel(displayResults,'Title','Various','Unit',...
    'normalized','Position',[0 .2 1 .19]);
HH0_handle = uicontrol(various_panel,'Style','checkbox','Unit','normalized',...
    'String','H/H_0','Position',[.05 .35 .3 .3]);
realPhi0 = uicontrol(various_panel,'Style','checkbox','Unit','normalized',...
    'String','real(phi_0)','Position',[.4 .35 .4 .3],'Value',1);

% Control buttons
uicontrol('Style','pushbutton','String','OK','Unit','normalized',...
    'Position',[0.25 .05 .375 .1],'Callback',@ok_Callback)
uicontrol('Style','pushbutton','String','Cancel','Unit','normalized',...
    'Position',[0.625 .05 .375 .1],'Callback',@cancel_Callback)


    %_________________
    function ok_Callback(hObject,eventdata)
    
    % Compute selected variables
    nOfTags = numel(tags);
    var2Plot = zeros(size(data.solution,1),nOfTags);
    title2Plot = cell(nOfTags,1);
    if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
        nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
        nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};
    elseif strcmp(data.computation,'NEFEM')
        nameNodal = 'X';
        nameCon = 'T';
    end
    
    k = 1;
    for itag = 1:nOfTags
        state = get(eval(tags{itag}),'Value');
        if state
            if itag == 1
                var2Plot(:,k) = real(data.solution);
                title2Plot{k} = 'real(\phi_r)';
            elseif itag == 2
                var2Plot(:,k) = imag(data.solution);
                title2Plot{k} = 'imag(\phi_r)';
            elseif itag == 3
                var2Plot(:,k) = abs(data.solution);
                title2Plot{k} = 'abs(\phi_r)';
            else
                if any(strcmp(data.computation,{'CDG' 'DG'}))
                    [aux,aux2] = totalPotential_CDG();
                else
                    aux2 = data.ip.value;
                    aux = aux2 + data.solution;
                end
                if itag == 4
                    var2Plot(:,k) = real(aux);
                    title2Plot{k} = 'real(\phi)';
                elseif itag == 5
                    var2Plot(:,k) = imag(aux);
                    title2Plot{k} = 'imag(\phi)';
                elseif itag == 6
                    var2Plot(:,k) = abs(aux);
                    title2Plot{k} = 'abs(\phi)';
                elseif itag == 7
                    var2Plot(:,k) = abs(aux)/data.ip.amplitude;
                    title2Plot{k} = 'H/H_0';
                elseif itag == 8
                    var2Plot(:,k) = real(aux2);
                    title2Plot{k} = 'real(\phi_0)';
                end
            end
            k = k + 1;
        end
    end
    
    % Plot selected variables
    for iplot = 1:k-1
        figure(2+iplot)
        if interpCheck
            degValueString = get(interpDeg_handle,'String');
            degValue = str2double(degValueString);
            plotHandle = plotSolutionInterpElem(data.mesh,var2Plot(:,iplot),...
                data.computation,degValue);
        else
            if strcmp(data.computation,'FEM')
                plotHandle = plotSolution(data.mesh.(nameNodal),data.mesh.(nameCon),...
                    var2Plot(:,iplot),data.mesh.referenceElement);
            elseif strcmp(data.computation,'NEFEM')
                plotHandle = plotSolution_NEFEM(data.mesh,var2Plot(:,iplot),-1);
            end
        end
        set(plotHandle(end),'location','EastOutside')
        title(title2Plot{iplot})
        axis off
    end
    
    guidata(handles.MainFigure,data);
    close(displayResults)
    
    end
    %_________________
    
    
    %_________________
    function cancel_Callback(hObject,eventdata)
    
    close(displayResults)
    
    end
    %_________________
    
    
    %_________________
    function interpCheck_Callback(hObject,eventdata)
    
    if get(hObject,'Value')
        interpCheck = true;
        set([interpDeg_handle interpDegText_handle],'Enable','on')
    else
        interpCheck = false;
        set([interpDeg_handle interpDegText_handle],'Enable','off')
    end
    
    end
    %_________________
    
    
    %_________________
    function [aux,aux2] = totalPotential_CDG()
    
    nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};
    nOfElementNodes = size(data.mesh.referenceElement.NodesCoord,1);
    nOfElements = size(data.mesh.(nameCon),1);
    aux = zeros(size(data.solution));
    aux2 = aux;
    for ielem = 1:nOfElements
        ind = (ielem-1)*nOfElementNodes+1:ielem*nOfElementNodes;
        aux2(ind) = data.ip.value(data.mesh.(nameCon)(ielem,:)).';
        aux(ind) = data.solution(ind) + aux2(ind);
    end
    
    end
    %_________________
    
end
        
        
        
        
        
        
        
        
        
        
        
        
        
        

