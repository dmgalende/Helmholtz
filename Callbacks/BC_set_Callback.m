function BC_set_Callback(hObject,eventdata,handles)

data = guidata(hObject);

currentCondition = data.currentBoundary.condition;
currentParam = data.BC.parameters{data.currentBoundary.value}{currentCondition};
selBoundaryName = data.mesh.boundaryNames{data.currentBoundary.value};
conditionNames = get(handles.BC_popup,'String');
selConditionName = conditionNames{currentCondition};

if any(currentCondition == [1 2 3])
    
    storeAndSet
    if strcmp(data.PML{5,data.currentBoundary.value},'on') %The boundary was a PML one
        set(handles.run_boundary,'Enable','off')
        set([handles.run_domain handles.run_all],'Value',0)
        data.PML{5,data.currentBoundary.value} = 'off'; %The boundary now is not a PML one
    end
    
elseif currentCondition == 4
    
    data.PML{5,data.currentBoundary.value} = 'on'; %The boundary now is a PML one
    if isempty(currentParam) || any(currentParam ~= data.currentBoundary.parameter)
        setOutput({'Computing PML absorption parameter...'},handles.run_wipOutput)
        
        %Compute absorption (sigma) parameter
        if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
            nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
        elseif strcmp(data.computation,'NEFEM')
            nameNodal = 'X';
        end
        param = data.currentBoundary.parameter;
        if param(4) > 1 && isempty(data.PML{6,data.currentBoundary.value})
            pos = param(4) - 1;
            nameSubCon = data.mesh.fieldNames{data.mesh.submeshIndex(pos)};
            PMLnodes = getNodes(data.mesh.(nameSubCon));
            data.PML{6,data.currentBoundary.value} = PMLnodes;
        elseif param(4) > 1 && ~isempty(data.PML{6,data.currentBoundary.value})
            PMLnodes = data.PML{6,data.currentBoundary.value};
        else
            PMLnodes = 1:size(data.mesh.(nameNodal),1);
            data.PML(6,data.currentBoundary.value) = {[]};
        end
        [sigma,sigmaPos,stretching] = coefPML(data.mesh.(nameNodal),PMLnodes,...
                                   data.PML{2,data.currentBoundary.value},...
                                   param(1),param(2),param(3));
        data.PML{1,data.currentBoundary.value} = zeros(size(data.mesh.(nameNodal),1),2);
        data.PML{1,data.currentBoundary.value}(PMLnodes,:) = sigma;
        data.PML{4,data.currentBoundary.value} = sigmaPos;
        data.PML{7,data.currentBoundary.value} = stretching;

        %Plot sigma value for each node
        delete([data.plotPMLnodes{:,data.currentBoundary.value}])
        data.plotPMLnodes(:,data.currentBoundary.value) = ...
            plotPMLnodes(data.mesh.(nameNodal)(PMLnodes,:),sigmaPos,handles.axesHandle);
    else
        setOutput({'The stored PML absorption parameter will be used'},handles.run_wipOutput)
    end
    storeAndSet
    set(handles.run_boundary,'Enable','off')
    set([handles.run_domain handles.run_all],'Value',0)
    set(handles.BC_PMLnodesCheck,'Enable','on','Value',0)
end

setOutput({['Condition on boundary "' selBoundaryName...
            '" set to "' selConditionName '"']},handles.run_wipOutput)

guidata(hObject,data);


    %_________________
    function storeAndSet
        
    data.BC.values(data.currentBoundary.value) = currentCondition;
    data.BC.parameters{data.currentBoundary.value}{currentCondition} =...
        data.currentBoundary.parameter;
    set(handles.run_boundary,'Value',0)
    set(handles.res_display,'Enable','off')
    if isfield(data,'solution')
        data = rmfield(data,'solution');
    end
    
    end
    %_________________

    
end


