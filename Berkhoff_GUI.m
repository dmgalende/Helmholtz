%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    USER INTERFACE - BERKHOFF EQUATION    %
% -----                                    %
% Developed by D. Modesto & G. Giorgiani   %
% January 2009                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clf, close all
home

%--------------------------------------------
% Main figure -  interface & general commands
%--------------------------------------------

MainFigure = figure('Visible','off','Name','Berkhoff Equation User Interface'...
    ,'MenuBar','none','Units','normalized','Position',[.15 .2 .6 .6],...
    'Tag','MainFigure');
axesHandle = axes('Position',[.03 .30 .67 .65],'Tag','axesHandle','Visible',...
    'off');

run setpath

%----------------------------
% Definition of the interface
%----------------------------

tags = {
        %----tools
        'zoomXY'                 2
        'panXY'                  2
        'fixBottom'              2
        'SimParam'               2
        'compOptions'            2
        %----computation
        'comp_panel'             3
        'comp_FEM'               0
        'comp_NEFEM'             0
        'comp_CDG'               0
        'comp_DG'                0
        %----mesh
        'mesh_panel'             0          
        'mesh_list'              1
        %----bottom
        'bottom_panel'           0          
        'bottom_list'            1
        %----boundary conditions
        'BC_panel'               0
        'BC_list'                1
        'BC_popup'               1
        'BC_set'                 1
        'BC_subpanel'            0
        'BC_paramText'           0
        'BC_value'               1
        'BC_PMLvalue_n'          1
        'BC_PMLvalue_R'          1
        'BC_PMLvalue_len'        1
        'BC_PMLvalue_area'       1
        'BC_PMLtext_n'           0
        'BC_PMLtext_R'           0
        'BC_PMLtext_len'         0
        'BC_PMLtext_area'        0
        'BC_PMLnodesCheck'       1  
        'BC_PMLareaCheck'        1
        %----incident potential
        'ip_panel'               0
        'ip_boundaryCheck'       1
        'ip_waveNumberBoundary'  1
        'ip_waveNumberValue'     1
        'ip_periodValue'         1
        'ip_directionValue'      1
        'ip_amplitudeValue'      1
        'ip_waveNumberText'      0
        'ip_periodText'          0
        'ip_directionText'       0
        'ip_amplitudeText'       0
        %----run
        'run_domain'             1
        'run_boundary'           1
        'run_all'                1
        'run_wipPanel'           0
        'run_wipOutput'          0
        %----results
        'res_display'            1
        'res_save'               1
        'res_load'               1
        %----selection info
        'sel_panel'              0
        'sel_meshText'           0
        'sel_bottomText'         0
        'sel_elemText'           0
        'sel_elem'               0
        'sel_mesh'               0
        'sel_bottom'             0
        'sel_meshCheck'          1
        'sel_bottomCheck'        1
        'sel_axesCheck'          1
        };        
    
%------------- Tools   
zoomXY = uitoggletool('TooltipString','Zoom','CData',...
    iconRead(fullfile('Callbacks/Icons','ZoomXY.gif')),'Separator','on',...
    'Enable','off');

panXY = uitoggletool('TooltipString','Pan','CData',...
    iconRead(fullfile('Callbacks/Icons','PanXY.gif')),'Separator','on',...
    'Enable','off');

fixBottom = uitoggletool('TooltipString','Fix bottom in PML','CData',...
    iconRead(fullfile('Callbacks/Icons','Fix.gif')),'Separator','on',...
    'Enable','off');

SimParam = uitoggletool('TooltipString','DG Simulation Parameters','CData',...
    iconRead(fullfile('Callbacks/Icons','Temp.gif')),'Separator','on',...
    'Enable','off');

compOptions = uitoggletool('TooltipString','Computational options','CData',...
    iconRead(fullfile('Callbacks/Icons','Options.png')),'Separator','on',...
    'Enable','off');
%--------------------------------------------------------------

%------------- Computation
comp_panel = uibuttongroup('Title','Select method','Unit','normalized',...
    'Position',[.72 .01 .27 .13]);

comp_FEM = uicontrol(comp_panel,'Style','radiobutton','Unit','normalized',...
    'String','FEM (Continuous Galerkin Finite Element Method)','Value',1,'Position',[.01 .8 .99 .2]);

comp_NEFEM = uicontrol(comp_panel,'Style','radiobutton','Unit','normalized',...
    'String','NEFEM (Nurbs Enhanced Finite Element Method)','Value',0,...
    'Position',[.01 .55 .99 .2]);

comp_CDG = uicontrol(comp_panel,'Style','radiobutton','Unit','normalized',...
    'String','CDG (Compact Discontinuos Galerkin FEM)','Value',0,...
    'Position',[.01 .3 .9 .2]);

comp_DG = uicontrol(comp_panel,'Style','radiobutton','Unit','normalized',...
    'String','DG (Discontinuos Galerkin Time Formulation)','Value',0,...
    'Position',[.01 .05 .9 .2]);
%--------------------------------------------------------------

%------------- Mesh list
mesh_panel = uipanel('Title','Open a mesh file',...
    'Position',[.72 .8 .13 .15]);

mesh_list = uicontrol(mesh_panel,'Style','listbox',...
    'String','no reading','Unit','normalized','Position',[0 0 1 1]);
%--------------------------------------------------------------

%------------- Bottom list
bottom_panel = uipanel('Title','Open a bottom file',...
    'Position',[.86 .8 .13 .15]);

bottom_list = uicontrol(bottom_panel,'Style','listbox',...
    'String','no reading','Unit','normalized','Position',[0 0 1 1],...
    'Enable','off');
%--------------------------------------------------------------

%------------- Boundary conditions settings
BC_panel = uipanel('Title','Select boundary conditions',...
    'Position',[.84 .25 .15 .5]);

BC_list = uicontrol(BC_panel,'Style','listbox','String',...
    'no mesh loaded','Unit','normalized','Position',[0 .58 1 .4],...
    'Enable','off');

BC_popup = uicontrol(BC_panel,'Style','popupmenu','String',...
    {'none','Partial absorption','Radiation','PML + radiation'},'Value',...
    1,'Unit','normalized','Position',[0 .57 1 .01],'Enable','off');

BC_set = uicontrol(BC_panel,'Style','pushbutton','String','Set BC',...
    'Unit','normalized','Position',[.0 .0 1 .1],'Enable','off');

BC_subpanel = uipanel('Parent',BC_panel,'Title','Parameters',...
    'Position',[0 .1 1 .42]);

BC_paramText = uicontrol(BC_subpanel,'Style','text','String','parameter',...
    'Unit','normalized','Position',[.45 .31 .55 .2]);

BC_value = uicontrol(BC_subpanel,'Style','edit','String','value',...
    'Unit','normalized','Position',[.05 .4 .4 .15],'Enable','off');

BC_PMLtext_n = uicontrol(BC_subpanel,'Style','text','String','Degree',...
    'Unit','normalized','Position',[.45 .72 .55 .2],'Visible','off');

BC_PMLvalue_n = uicontrol(BC_subpanel,'Style','edit','String','value',...
    'Unit','normalized','Position',[.05 .83 .4 .12],'Visible','off');

BC_PMLtext_R = uicontrol(BC_subpanel,'Style','text','String','R/Length',...
    'Unit','normalized','Position',[.45 .55 .55 .2],'Visible','off');

BC_PMLvalue_R = uicontrol(BC_subpanel,'Style','edit','String','value',...
    'Unit','normalized','Position',[.05 .665 .4 .12],'Visible','off');

BC_PMLtext_len = uicontrol(BC_subpanel,'Style','text','String','Length',...
    'Unit','normalized','Position',[.45 .39 .55 .2],'Visible','off');

BC_PMLvalue_len = uicontrol(BC_subpanel,'Style','edit','String','value',...
    'Unit','normalized','Position',[.05 .5 .4 .12],'Visible','off');

BC_PMLtext_area = uicontrol(BC_subpanel,'Style','text','String',...
    'Subdomain','Unit','normalized','Position',[.45 .22 .55 .2],...
    'Visible','off');

BC_PMLvalue_area = uicontrol(BC_subpanel,'Style','popupmenu','String',...
    'none','Unit','normalized','Position',[.05 .25 .4 .2],'Value',1,...
    'Visible','off');

BC_PMLnodesCheck = uicontrol(BC_subpanel,'Style','checkbox','String',...
    'View PML nodes','Unit','normalized','Position',[.05 .02 .65 .1],...
    'Visible','off','Enable','off');

BC_PMLareaCheck = uicontrol(BC_subpanel,'Style','checkbox','String',...
    'View PML area','Unit','normalized','Position',[.05 .15 .65 .1],...
    'Visible','off','Enable','off');
%--------------------------------------------------------------

%------------- Incident potential settings
ip_panel = uipanel('Title','Enter incident potential data',...
    'Position',[.03 .01 .19 .25]);

ip_boundaryCheck = uicontrol(ip_panel,'Style','checkbox','Unit',...
    'normalized','Position',[.05 .84 .95 .1],'String',...
    'Wave number from a boundary','Enable','off');

ip_waveNumberBoundary = uicontrol(ip_panel,'Style','popupmenu',...
    'Unit','normalized','String','no reading','Position',[.05 .64 .4 .1],...
    'Visible','off');

ip_waveNumberValue = uicontrol(ip_panel,'Style','edit',...
    'Unit','normalized','String','value','Position',[.05 .64 .4 .1]);

ip_periodValue = uicontrol(ip_panel,'Style','edit','String',...
    'value','Unit','normalized','Position',[.05 .45 .4 .1]);

ip_directionValue = uicontrol(ip_panel,'Style','edit','String',...
    'value','Unit','normalized','Position',[.05 .26 .4 .1]);

ip_amplitudeValue = uicontrol(ip_panel,'Style','edit','String',...
    'value','Unit','normalized','Position',[.05 .07 .4 .1]);

ip_waveNumberText = uicontrol(ip_panel,'Style','text','Unit','normalized',...
    'Position',[.5 .52 .45 .2],'String','wave number');

ip_periodText = uicontrol(ip_panel,'Style','text','Unit','normalized',...
    'Position',[.5 .33 .45 .2],'String','period [sec]');

ip_directionText = uicontrol(ip_panel,'Style','text','Unit','normalized',...
    'Position',[.5 .24 .45 .1],'String','direction [deg]');

ip_amplitudeText = uicontrol(ip_panel,'Style','text','Unit','normalized',...
    'Position',[.5 .055 .45 .1],'String','amplitude');
%--------------------------------------------------------------

%------------- Run
run_domain = uicontrol('Style','toggleButton','String','Run domain',...
    'Unit','normalized','Position',[.28 .2 .1 .05],'Enable','off');

run_boundary = uicontrol('Style','toggleButton','String','Run boundary',...
    'Unit','normalized','Position',[.28 .15 .1 .05],'Enable','off');

run_all = uicontrol('Style','toggleButton','String','Run all',...
    'Unit','normalized','Position',[.44 .2 .1 .05],'Enable','off');

run_wipPanel = uipanel('Title','Output','Position',[.28 .01 .42 .13]);

run_wipOutput = uicontrol(run_wipPanel,'Style','listbox','String',...
    {'Ready'},'Unit','normalized','Position',[0 0 1 1]);
%--------------------------------------------------------------

%------------- Results
res_display = uicontrol('Style','pushbutton','String','Display results',...
    'Unit','normalized','Position',[.44 .15 .1 .05],'Enable','off');

res_save = uicontrol('Style','pushbutton','String','Save',...
    'Unit','normalized','Position',[.6 .2 .1 .05],'Enable','off');

res_load = uicontrol('Style','pushbutton','String','Load',...
    'Unit','normalized','Position',[.6 .15 .1 .05]);
%--------------------------------------------------------------

%------------- Selection info
sel_panel = uipanel('Title','Selection Info','Position',[.72 .4 .11 .35]);

sel_meshText = uicontrol(sel_panel,'Style','text','String',...
    'Mesh file:','Unit','normalized','HorizontalAlignment','left',...
    'Position',[0 .78 1 .2]);

sel_bottomText = uicontrol(sel_panel,'Style','text','String',...
    'Bottom file:','Unit','normalized','HorizontalAlignment','left',...
    'Position',[0 .58 1 .2]);

sel_mesh = uicontrol(sel_panel,'Style','text','String',...
    'none','Unit','normalized','FontWeight','bold','Position',...
    [0 .8 1 .1]);

sel_bottom = uicontrol(sel_panel,'Style','text','String',...
    'none','Unit','normalized','FontWeight','bold','Position',...
    [0 .6 1 .1]);

sel_meshCheck = uicontrol(sel_panel,'Style','checkbox','Unit',...
    'normalized','Position',[.05 .14 .7 .1],'String',...
    'View mesh','Enable','off','Value',0);

sel_bottomCheck = uicontrol(sel_panel,'Style','checkbox','Unit',...
    'normalized','Position',[.05 .01 .73 .1],'String',...
    'View bottom','Enable','off','Value',0);

sel_axesCheck = uicontrol(sel_panel,'Style','checkbox','Unit',...
    'normalized','Position',[.05 .27 .67 .1],'String',...
    'View axes','Value',0);

sel_elemText = uicontrol(sel_panel,'Style','text','String',...
    'Element Type:','Unit','normalized','HorizontalAlignment','left',...
    'Position',[0 .38 1 .2]);

sel_elem = uicontrol(sel_panel,'Style','text','String',...
    'none','Unit','normalized','FontWeight','bold','Position',...
    [0 .4 1 .1]);
%--------------------------------------------------------------

%--------------------------------------
% Initialization, variables & callbacks
%--------------------------------------

% Assign tags
nTags = size(tags,1);
for item = 1:nTags
    itagString = tags{item,1};
    itag = eval(itagString);
    set(itag,'Tag',itagString)
end

% Variables & Default Settings
dir_struct = struct('mesh',struct,'bottom',struct);
ip = struct('waveNumberValue',[],'period',[],'direction',[],'amplitude',[]);
sp = struct('sstBoundary',1,'sstol',1e-3,'nOfCyclesMax',30,'stabParam',4,...
    'spaceXtimeFrame',0.04,'optStore',0);
cpuTimeInfo = struct('volMat',[],'volMatToc',[],'boundary',[],'boundaryToc',[],...
    'linearSystem',[],'linearSystemToc',[]);
data = struct('dir_struct',dir_struct,'ip',ip,'sp',sp,'cputime',cpuTimeInfo,...
    'computation','FEM','staticCondensation',false,'constantBottomFlag',false);
guidata(MainFigure,data);
handles = guihandles;

% Assign callbacks
for item = 1:nTags
    itagString = tags{item,1};
    icall = tags{item,2};
    if icall
        itag = eval(itagString);
        ifun = str2func([itagString '_Callback']);
        if icall == 1
            set(itag,'Callback',{ifun,handles})
        elseif icall == 2
            set(itag,'ClickedCallback',{ifun,handles})
        elseif icall == 3
            set(itag,'SelectionChangeFcn',{ifun,handles})
        end
    end
end

% Initialization tasks
mesh_dir = [pwd '/Meshes'];
bottom_dir = [pwd '/Bottoms'];
load_listbox({mesh_dir bottom_dir},...
    [handles.mesh_list handles.bottom_list],{'mesh' 'bottom'})
set(MainFigure,'Visible','on')

