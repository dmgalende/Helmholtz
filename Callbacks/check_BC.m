function check_BC(value,handles)

data = guidata(handles.MainFigure);

param = data.BC.parameters{data.currentBoundary.value}{value};
handlesPML = [handles.BC_PMLvalue_n
              handles.BC_PMLvalue_R
              handles.BC_PMLvalue_len
              handles.BC_PMLvalue_area
              handles.BC_PMLareaCheck
              handles.BC_PMLnodesCheck
              handles.BC_PMLtext_area
              handles.BC_PMLtext_len
              handles.BC_PMLtext_n
              handles.BC_PMLtext_R];
switch value
    case 1
        data.currentBoundary.condition = value;
        set(handles.BC_paramText,'String','parameter','Visible','on')
        set(handles.BC_value,'String','value','Enable','off','Visible','on')
        set(handles.BC_set,'Enable','on')
        set(handlesPML,'Visible','off')
        data.currentBoundary.parameter = [];
    case 2
        data.currentBoundary.condition = value;
        set(handles.BC_paramText,'String','alpha','Visible','on')
        set(handles.BC_value,'Enable','on','Visible','on')
        
        if isempty(param) 
            set(handles.BC_set,'Enable','off')
            set(handles.BC_value,'String','value')
            data.currentBoundary.parameter = [];
        else
            set(handles.BC_set,'Enable','on')
            set(handles.BC_value,'String',num2str(param))
            data.currentBoundary.parameter = param;
        end
        
        set(handlesPML,'Visible','off')
    case 3
        data.currentBoundary.condition = value;
        set(handles.BC_paramText,'String','radius','Visible','on')
        set(handles.BC_value,'Enable','on','Visible','on')
        
        if isempty(param) 
            set(handles.BC_set,'Enable','off')
          	set(handles.BC_value,'String','value')
            data.currentBoundary.parameter = [];
        else
            set(handles.BC_set,'Enable','on')
            set(handles.BC_value,'String',num2str(param))
            data.currentBoundary.parameter = param;
        end
        
        set(handlesPML,'Visible','off')
    case 4
        
        %First, check the posibility of assigning the PML condition to the
        %current boundary and return to the previous condition in negative
        %case
        if isempty(data.PML{2,data.currentBoundary.value})
            try
                nameBoundary = data.mesh.fieldNames{data.mesh.boundaryIndex(data.currentBoundary.value)};
                if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
                    nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
                    [data.PML{2,data.currentBoundary.value} data.PML{3,data.currentBoundary.value}] = ...
                    getPMLboundaries(data.mesh.(nameNodal),data.mesh.(nameBoundary),...
                                     data.mesh.referenceElement,...
                                     handles.run_wipOutput);
                elseif strcmp(data.computation,'NEFEM')
                    [data.PML{2,data.currentBoundary.value} data.PML{3,data.currentBoundary.value}] = ...
                    getPMLboundaries(data.mesh.X,data.mesh.(nameBoundary),...
                                     data.mesh.referenceElement,...
                                     handles.run_wipOutput,...
                                     data.mesh.nurbs,...
                                     data.mesh.trimmedInfo.(nameBoundary));
                end
                data.currentBoundary.condition = value;
            catch %Execute these lines if the boundary is not a cartesian one
                data.currentBoundary.condition = data.BC.values(data.currentBoundary.value);
                set(handles.BC_popup,'Value',data.currentBoundary.condition)
                check_BC(data.currentBoundary.condition,handles)
                error('Error in Berkhoff GUI')
            end
        else
            data.currentBoundary.condition = value;
        end                  
        
        set([handles.BC_paramText handles.BC_value],'Visible','off')
        set(handlesPML,'Visible','on')
        data.currentBoundary.parameter = [-1,-1,-1,1]; 
                                        
        if ~isempty(param)
            for cont = 1:3
                set(handlesPML(cont),'String',num2str(param(cont)))
                data.currentBoundary.parameter(cont) = param(cont);
            end
            set(handlesPML(4),'Value',param(4))
            set([handles.BC_set handlesPML(5)],'Enable','on')
            if param(3) ~= data.plotPMLarea{2,data.currentBoundary.value}
                data.plotPMLarea{2,data.currentBoundary.value} = param(3);
                delete(data.plotPMLarea{1,data.currentBoundary.value})          
                data.plotPMLarea{1,data.currentBoundary.value} = ...
                    plotPMLarea(data.PML{2,data.currentBoundary.value},param(3),handles.axesHandle);
            end
            data.currentBoundary.parameter(4) = param(4);
        else
            for cont = 1:3
                set(handlesPML(cont),'String','value')
            end
            set(handlesPML(4),'Value',1)
            set([handles.BC_set handlesPML(5)],'Enable','off')
        end
        
        set(handlesPML(6),'Enable',data.PML{5,data.currentBoundary.value})
end

checkVisibility

guidata(handles.MainFigure,data);


    %_________________
    function checkVisibility
    
    %Submeshes
    if numel(data.mesh.submeshNames) > 0
        set(data.plotSubmesh,'Visible','off')
        if value == 4 && ~isempty(param) && param(4) > 1
            set(data.plotSubmesh(param(4)-1),'Visible','on')
        end
    end
    
    %PML area and nodes (always not visible)
    for aux = 1:numel(data.mesh.boundaryNames)
        set(data.plotPMLarea{1,aux},'Visible','off')
        set([data.plotPMLnodes{:,aux}],'Visible','off')
    end
    set(handlesPML([5 6]),'Value',0)
    
    end
    %_________________
    
    
end
            
                
        
        
        
        
        
        
        
        
        
        
