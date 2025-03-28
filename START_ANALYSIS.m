clc
clear all

% put back cursors to figure 
set(groot, ...
    'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot, ...
    'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

% create main window, search if the figure exists
hf = findall(0,'tag','selectAnalysis');
if isempty(hf)
    % Launch the figure
    % main_figure
    mainFig = figure('Name','select analysis','units','normalized', ...
        'Position',[0.375 0.25 0.25 0.5]);
    set(mainFig, 'PaperPositionMode','auto',...
        'PaperOrientation','landscape',...
        'PaperType','A4',...
        'Tag','selectAnalysis');
else
    % Figure exists so bring Figure to the focus
    figure(hf);
end

% remove old paths
try  
    pathStr = path;
    C = strsplit(pathStr,':')';
    m_path = cellfun(@(x) isempty(x), regexp(C,'/MATLAB_R20') );
    pathStrRm = strjoin(C(m_path),':');
    rmpath(pathStrRm);
catch
    pathStr = path;
    C = strsplit(pathStr,':')';
    m_path = cellfun(@(x) isempty(x), regexp(C,'\\MATLAB\\R20') );
    pathStrRm = strjoin(C(m_path),':');
    rmpath(pathStrRm);
end

% set working directory
mfilePath = mfilename('fullpath');
wd_path = fileparts(mfilePath);
addpath(wd_path) % do not include subfolders
%addpath(genpath(currentFolder))
% add path of folders containing functions for loading image
addpath(fullfile(wd_path,'OME_bioformats'))
addpath(fullfile(wd_path,'sharedFunctions'))
addpath(fullfile(wd_path,'Igor2Matlab'))

if exist("mainFig", "var")
    % panel for choose analysis
    hp1 = uipanel('Title','select type of calcium signals:', ...
        'Parent',mainFig, 'Position',[0.05 0.05 0.9 0.9],...
        'Visible','on', 'Units','normalized', 'FontUnits','normalized',...
        'FontSize',0.05, 'FontWeight','bold');

    w = 0.9;
    h = 0.25;
    N_pb = 3;
    nR = 3;
    nC = 1;
    dy = (1-nR*h)/(nR+1);
    dx = (1-nC*w)/(nC+1);

    pb_names = {'calcium sparks detection', ...
        'calcium sparks recovery',...
        'calcium transients & waves'};
    c = {'k','b','r'};

    for i=1:nR*nC

        col = ceil(i/nR);
        row = mod(i,nR); if row==0, row = nR; end

        uicontrol('Style','pushbutton','String',pb_names{i},...
            'FontUnits','normalized','Parent',hp1,'Units','normalized',...
            'FontSize',0.225,'Position',[col*dx+(col-1)*w 1-row*h-row*dy w h],...
            'FontWeight','bold','ForegroundColor',c{i},...
            'Callback',{@selectAnalysis, wd_path});

        if i==N_pb
            break
        end

    end
end



% h_push_SpDet = uicontrol('Style','pushbutton','String','spark detection',...
%     'FontUnits','normalized','Parent',hp1,'Units','normalized',...
%     'FontSize',0.25,'Position',[dx (nR-1)*h+(nR-0)*dy w h],...
%     'FontWeight','bold','ForegroundColor','k',...
%     'Callback',@selectAnalysis);
% 
% h_push_SpRecRy = uicontrol('Style','pushbutton','String','spark recovery ryanodine',...
%     'FontUnits','normalized','Parent',hp1,'Units','normalized',...
%     'FontSize',0.25,'Position',[dx (nR-2)*h+(nR-1)*dy w h],...
%     'FontWeight','bold','ForegroundColor','b',...
%     'Callback',@selectAnalysis);
% 
% h_push_SpRecP = uicontrol('Style','pushbutton','String','spark recovery photolysis',...
%     'FontUnits','normalized','Parent',hp1,'Units','normalized',...
%     'FontSize',0.25,'FontWeight','bold','ForegroundColor','r',...
%     'Position',[dx (nR-3)*h+(nR-2)*dy w h],...
%     'Callback',@selectAnalysis);
% 
% h_push_TransWave = uicontrol('Style','pushbutton','String','trasients & waves',...
%     'FontUnits','normalized','Parent',hp1,'Units','normalized',...
%     'FontSize',0.25,'FontWeight','bold','ForegroundColor','c',...
%     'Position',[nC*dx+(nC-1)*w (nR-1)*h+(nR-0)*dy w h],...
%     'Callback',@selectAnalysis);
% 
% h_push_SparkBlinkPairs = uicontrol('Style','pushbutton','String','spark-blink pairs',...
%     'FontUnits','normalized','Parent',hp1,'Units','normalized',...
%     'FontSize',0.25,'FontWeight','bold','ForegroundColor','m',...
%     'Position',[nC*dx+(nC-1)*w (nR-2)*h+(nR-1)*dy w h],...
%     'Callback',@selectAnalysis);





