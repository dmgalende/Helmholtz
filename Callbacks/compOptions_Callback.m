function compOptions_Callback(hObject,eventdata,handles)

data = guidata(handles.MainFigure);

% Options figure
displayOptions = figure('Name','Computational options','MenuBar','none','Unit',...
    'normalized','Position',[.35 .5 .3 .2],'WindowStyle','modal');

% Options
panel1 = uipanel(displayOptions,'Unit','normalized','Position',[0 .175 1 .825]);
static = uicontrol(panel1,'Style','checkbox','Unit','normalized',...
    'String','Static condensation','Position',[.05 .8 .4 .1],'Value',data.staticCondensation);
constantBottom = uicontrol(panel1,'Style','checkbox','Unit','normalized',...
    'String','Constant bottom','Position',[.05 .6 .4 .1],'Value',data.constantBottomFlag);

% Control buttons
uicontrol('Style','pushbutton','String','OK','Unit','normalized',...
    'Position',[0.25 .05 .375 .1],'Callback',@ok_Callback)
uicontrol('Style','pushbutton','String','Cancel','Unit','normalized',...
    'Position',[0.625 .05 .375 .1],'Callback',@cancel_Callback)


    %_________________
    function ok_Callback(hObject,eventdata)
        
        %Static condensation
        if get(static,'value') && data.mesh.referenceElement.degree < 3
            setOutput({'Element degree is less than cubic, static condensation wont be done'},...
                handles.run_wipOutput)
            set(static,'value',0)
        end
        data.staticCondensation = get(static,'value');
        %Constant bottom
        data.constantBottomFlag = get(constantBottom,'value');
        
        guidata(handles.MainFigure,data)
        set(handles.compOptions,'state','off')
        close(displayOptions)
        
    end
    %_________________
    
    %_________________
    function cancel_Callback(hObject,eventdata)
        
        set(handles.compOptions,'state','off')
        close(displayOptions)
        
    end
    %_________________
   
    
end
        