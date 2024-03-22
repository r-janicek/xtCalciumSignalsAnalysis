%% create main figure and axes   
mainFig = figure('Name','transients and waves analysis', ...
    'units','normalized', 'outerposition',[0 0.05 1 0.95]);
set(mainFig, 'PaperPositionMode','auto', ...
             'PaperOrientation','landscape', ...
             'PaperType','A4', ...
             'Tag','mainFig');
% create axes
ax_img = axes('Parent',mainFig);
set(ax_img,'Position',[0.03 0.57 0.94 0.4])
% set(get(ax_img,'Xlabel'),'String','t (ms)')
set(ax_img,'XTick',[],'FontSize',14)
set(get(ax_img,'Ylabel'),'String','x (pixels)','FontWeight','bold')

% set up scale
ax_sc = axes('Parent',mainFig); % invisible axes
set(ax_sc, 'Visible', 'off','Position',[0.975 0.7075 0.005 0.125]);
% Switch off autoscaling.
set(ax_sc, 'Xlim', [0, 1], 'YLim', [0, 1]);
line([0.5 0.5],[0,1],'Parent',ax_sc,'LineWidth',4,'Color','k')
h_txt_scale_img = text('Parent',ax_sc, 'Position',[1 0.5], ...
    'String',[sprintf('%0.2f',0),' \mum'],...
    'FontSize',16, 'FontWeight','bold', 'Rotation',90,...
    'HorizontalAlignment','center', 'VerticalAlignment','cap');

ax_prof = axes('Parent',mainFig,'units','normalized','PickableParts','all');
set(ax_prof,'Position',[0.03 0.36 0.94 0.2],'YGrid','on','FontSize',14)

set(get(ax_prof,'Xlabel'),'String','t (ms)','FontWeight','bold')
set(get(ax_prof,'Ylabel'),'String','Fluorescence (F)','FontWeight','bold')
% mouse baseline selection tool
set(ax_prof,'buttondownfcn',@mouseSetMaskFcn)

hlink1 = linkprop([ax_img,ax_prof],'XLim');
%hlink2 = linkprop([ax_img_sparks,ax_img_sparks_2],'YLim');
hlink2 = [];

%linkaxes([ax_img,ax_img_sparks,ax_img_orig],'xy')
%linkaxes([ax_img,ax_img_sparks,ax_prof,ax_img_orig],'x')

hObjs.mainFig = mainFig;
hObjs.ax_img = ax_img;
hObjs.ax_sc = ax_sc;
hObjs.h_txt_scale_img = h_txt_scale_img;
hObjs.ax_prof = ax_prof;
hObjs.hlink1 = hlink1;
hObjs.hlink2 = hlink2;


%% text to show image path and name and channels assignment 
txt_name_img = uicontrol('Style','text',...
    'Parent',mainFig,'Units','normalized','Position', [0.25 0.980 0.51 0.02],...
    'FontUnits','normalized','FontSize',0.65,'FontWeight','bold',...
    'HorizontalAlignment','left','String','image name');

% individual channels assignment 
uicontrol('Style','text',...
    'Parent',mainFig,'Units','normalized','Position', [0.01 0.980 0.03 0.02],...
    'FontUnits','normalized','FontSize',0.75,'FontWeight','bold',...
    'HorizontalAlignment','left','String','fluo ch:');
popUpMenuChFluo = uicontrol('Style','popup',...
    'Parent',mainFig,'Units','normalized','Position', [0.04 0.980 0.05 0.02],...
    'FontUnits','normalized','FontSize',0.75,'Tag','mainImgFluoCh',...
    'String',{'ch #1','ch #2','none'},'Callback', {@setImgChannel,mainFig});

uicontrol('Style','text',...
    'Parent',mainFig,'Units','normalized','Position', [0.1 0.980 0.04 0.02],...
    'FontUnits','normalized','FontSize',0.75,'FontWeight','bold',...
    'HorizontalAlignment','left','String','trans ch:');
popUpMenuChTrans = uicontrol('Style','popup',...
    'Parent',mainFig,'Units','normalized','Position', [0.14 0.980 0.05 0.02],...
    'FontUnits','normalized','FontSize',0.75,'Tag','mainImgTransCh',...
    'String',{'ch #1','ch #2','none'},'Callback', {@setImgChannel,mainFig});

% axes to show trans channel
h_ax_transCh = axes('Parent',mainFig);
set(h_ax_transCh,'Position',[0.19 0.975 0.05 0.0209])
set(h_ax_transCh,'XTick',[],'YTick',[],'FontSize',8)

hObjs.txt_name_img = txt_name_img;
hObjs.h_ax_transCh = h_ax_transCh;
hObjs.popUpMenuChFluo = popUpMenuChFluo;
hObjs.popUpMenuChTrans = popUpMenuChTrans;


%% create slider when image is too big in time direction

h_bg_sld = uibuttongroup('Visible','on',...
    'Parent',mainFig,'Units','normalized','Position', [0.8 0.975 0.125 0.025],...
                  'SelectionChangedFcn',{@bg_sld_selection});
                     
% Create three radio buttons in the button group.
h_b1 = uicontrol(h_bg_sld,'Style','radiobutton','Units','normalized',...
                  'String','whole image','FontUnits','normalized','FontSize',0.75,...
                  'Position',[0.05 0.05 0.55 0.9],...
                  'HandleVisibility','off');
              
h_b2 = uicontrol(h_bg_sld,'Style','radiobutton','Units','normalized',...
                  'String','#s part','FontUnits','normalized','FontSize',0.75,...
                  'Position',[0.6 0.05 0.45 0.9],...
                  'HandleVisibility','off');

uicontrol('Style','text','Parent',mainFig,'Units','normalized',...
            'Position', [0.96 0.975 0.01 0.025],...
            'FontUnits','normalized','FontSize',0.75,'String','s');
h_edit_imgSld_w = uicontrol('Style','edit','String','5',...
                    'FontUnits','normalized','FontSize',0.75,...
                    'Parent',mainFig,'Units','normalized',...
                    'Position',[0.935 0.975 0.02 0.025],...
                    'Callback', {@bg_sld_selection});
                
h_img_sld = uicontrol('Style', 'slider','Parent',mainFig,'FontUnits','normalized',...
                'Min',1,'Max',100,'Value',1,'Units','normalized',...
                'Position', [0.03 0.56 0.94 0.01],'BackgroundColor',[0.7 0.7 0.7],...
                'Callback', '','Enable','off'); %{@imgShow_slider}

       
hObjs.h_bg_sld = h_bg_sld;
hObjs.h_b1_sld = h_b1;
hObjs.h_b2_sld = h_b2;
hObjs.h_edit_imgSld_w = h_edit_imgSld_w;
hObjs.h_img_sld = h_img_sld;


%% open image and set ROI for cropping
h_pb_open = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> open <br> image <html>',...
    'FontUnits','normalized','FontSize',0.225,'FontWeight','bold',...
    'Parent',mainFig,'Units','normalized','Position', [0.01 0.2 0.05 0.1],...
    'Callback',{@openImg,mainFig},'Enable','on');

% crop image
h_pb_crop = uicontrol('Style','pushbutton',...
    'String','<html> <p align="center"> set ROI <br> for cropping <html>',...
    'FontUnits','normalized', 'FontSize',0.2,...
    'Parent',mainFig, 'Units','normalized', ...
    'Position', [0.01 0.09 0.05 0.1],...
    'Callback',{@cropImg,mainFig}, 'Enable','on');

% % baseline correction
% h_pb_bsCorr = uicontrol('Style', 'pushbutton',...
%     'String','baseline correction','FontUnits','normalized','FontSize',0.4,...
%     'Parent',mainFig,'Units','normalized','Position', [0.03 0.01 0.1 0.05],...
%     'Callback',{@baselineCorrection,mainFig},'Enable','on');

hObjs.h_pb_open = h_pb_open;
hObjs.h_pb_crop = h_pb_crop;


%% pixel sizes panel
hp_pxSz = uipanel('Title','pixel sizes:', 'Parent',mainFig,...
    'Position',[0.01 0.02 0.1 0.06], 'FontUnits','normalized',...
    'FontSize',0.25, 'FontWeight','bold');

uicontrol('Style','text', 'Parent',hp_pxSz, 'Units','normalized',...
        'Position',[0.01 0.05 0.2 0.9], ...
        'String',{'x', '(µm)'}, ...
        'FontUnits','normalized','FontSize',0.4);       
h_edit_pxSzX = uicontrol('Style','edit', 'String','1',...
        'Parent',hp_pxSz,'Units','normalized', ...
        'Position',[0.22 0.05 0.255 0.9], 'Tag','pxSzX',...
        'FontUnits','normalized','FontSize',0.4, ...
        'Callback', {@setPxSzManually, mainFig});
uicontrol('Style','text', 'Parent',hp_pxSz, 'Units','normalized',...
        'Position',[0.525 0.05 0.2 0.9], ...
        'String',{'t', '(ms)'}, ...
        'FontUnits','normalized','FontSize',0.4);       
h_edit_pxSzT = uicontrol('Style','edit', 'String','1',...
        'Parent',hp_pxSz,'Units','normalized', ...
        'Position',[0.735 0.05 0.255 0.9], 'Tag','pxSzT',...
        'FontUnits','normalized','FontSize',0.4, ...
        'Callback', {@setPxSzManually, mainFig});
hObjs.h_edit_pxSzX = h_edit_pxSzX;
hObjs.h_edit_pxSzT = h_edit_pxSzT;


%% animal/experiment notes
hpA = uipanel('Title','animal & notes','Parent',mainFig,'Position',[0.07 0.09 0.125 0.21],...
    'Visible', 'on','Units','normalized','FontUnits','normalized',...
    'FontSize',0.1,'FontWeight','bold');

popUpMenuAnimal = uicontrol('Style','popup',...
    'Parent',hpA,'Units','normalized','Position', [0.01 0.9 0.98 0.1],...
    'FontUnits','normalized','FontSize',0.8,'Tag','popUpMenuAnimal',...
    'String',{'wt',...
              'PA-RFP RyR', ...
              'wt (S2030A littermate)',...
              'wt (R420Q littermate)',...
              'S2030A',...
              'S2808A/14A',...                         
              'R420Q (homozygous)',...
              'R420Q (heterozygous)', ...
              'IP3R/tTA', ...
              'wt/tTA', ...
              'unknown'}); 

ScSz = get(0,'ScreenSize');

h_table_notes = uitable( ...
    'Data',{'ctrl'; ...
            'animal #: '; ...
            '1.8 mM CaCl2 tyrode'; ...
            'field stim.: 30 s, 1 Hz (0.5 ms, 30 V)'; ...
            'loading: 2 uM Cal520 AM, 60 min (RT)';''},...
    'ColumnEditable',true,...
    'ColumnName',{'<HTML> <font size="4"> <b> notes </b> </HTML>'},...
    'ColumnWidth',{hpA.Position(3)*0.8*ScSz(3)},'Parent',hpA,...
    'Units','normalized','Position',[0.01 0.05 0.98 0.8],...
    'FontUnits','normalized','FontSize',0.1);

hObjs.popUpMenuAnimal = popUpMenuAnimal;
hObjs.h_table_notes = h_table_notes;


%% panel for blank
hp3 = uipanel('Title','blank','Parent',mainFig,'Position',[0.205 0.09 0.075 0.21],...
    'Visible', 'on','Units','normalized','FontUnits','normalized',...
    'FontSize',0.1,'FontWeight','bold');

h_pb_GetBlank = uicontrol('Style', 'pushbutton',...
    'String','cell (xy img)','FontUnits','normalized','FontSize',0.3,...
    'Parent',hp3,'Units','normalized','Position', [0.15 0.64 0.7 0.3],...
    'Callback', {@getBlankXYimg,mainFig},'Enable','on');

h_push_BlankROI = uicontrol('Style', 'pushbutton','String','set ROI for blank',...
    'FontUnits','normalized','Parent',hp3,'Units','normalized','FontSize',0.3,...
    'Position', [0.05 0.28 0.9 0.3],'Callback', {@setBlankROI,mainFig});

h_edit_Blank = uicontrol('Style', 'edit','String','nan',...
    'FontUnits','normalized','Parent',hp3,'Units','normalized',...
    'FontSize',0.5,'Position', [0.15 0.02 0.7 0.2],...
    'Callback', {@calculateBlank,mainFig});
    
hObjs.h_edit_Blank = h_edit_Blank;
hObjs.h_push_BlankROI = h_push_BlankROI;
hObjs.h_pb_GetBlank = h_pb_GetBlank;


%% panel normalization 
hp5 = uipanel('Title','self-ratio','Parent',mainFig,'Position',[0.29 0.09 0.1 0.21],...
    'FontUnits','normalized','FontSize',0.1,'FontWeight','bold');

%select baseline function
check_doBsFit = uicontrol('Style', 'checkbox','Parent',hp5,...
    'FontUnits','normalized','Value',0,'Units','normalized','FontSize',0.9,...
    'String','do baseline fitting','Position', [0.05 0.875 0.95 0.1],...
    'Tag','check_bsFit','Callback',{@checkboxFcn,mainFig});

popUpMenuBaselineFcn = uicontrol('Style','popup',...
    'Parent',hp5,'Units','normalized','Position', [0.05 0.75 0.90 0.1],...
    'FontUnits','normalized','FontSize',0.9,'Tag','popUpMenuBsFcn',...
    'String',{'exp1', 'exp2', 'stretchedExp', 'poly1', 'poly2',...
    'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'poly9'},...
    'Callback',{@selectBaselineFcn,mainFig},'Value',7,'Enable','off');

% add slider to control baseline percentile 
startVal = 20;
h_txt_sld_Bs = uicontrol('Style','text',...
    'Parent',hp5,'Units','normalized','Position', [0.01 0.6 0.98 0.1],...
    'FontUnits','normalized','FontSize',0.5,'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'String',sprintf('mask of baseline as %dth percentile',startVal));

% set up slider, percentiles 
sldMin = 1;
sldMax = 99; 

sldStep = 1/(sldMax-sldMin);
sld_Bs = uicontrol('Style', 'slider','Parent',hp5,'FontUnits','normalized',...
    'Min',sldMin,'Max',sldMax,'Value',startVal,'SliderStep',[sldStep 10*sldStep],...
    'Units','normalized','Position', [0.05 0.5 0.9 0.1],...
    'Callback', {@bsMask_slider,mainFig},'Enable','off');

% do normalization
h_pb_norm = uicontrol('Style', 'pushbutton','FontUnits','normalized',...
    'String','normalize data','FontSize',0.5,...
    'Parent',hp5,'Units','normalized','Position', [0.05 0.25 0.9 0.2],...
    'Callback', {@norm_data,mainFig});

h_pb_undo_norm = uicontrol('Style', 'pushbutton','String','start over',...
    'FontUnits','normalized','FontSize',0.5,...
    'Parent',hp5,'Units','normalized','Position', [0.05 0.01 0.9 0.2],...
    'Callback', {@undo_norm,mainFig});

hObjs.h_pb_norm = h_pb_norm;
hObjs.h_pb_undo_norm = h_pb_undo_norm;
hObjs.popUpMenuBaselineFcn = popUpMenuBaselineFcn;
hObjs.check_doBsFit = check_doBsFit;
hObjs.h_txt_sld_Bs = h_txt_sld_Bs;
hObjs.sld_Bs = sld_Bs;


%% panel for event image selection 
hp_SD = uipanel('Title','select & analyze events', ...
    'Parent',mainFig,...
    'Position',[0.4 0.09 0.34 0.21], ...
    'FontUnits','normalized',...
    'FontSize',0.1,'FontWeight','bold');

h_push_eventImgROI = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> ROI to select image <br> of event <html>',...
    'FontUnits','normalized', 'FontSize',0.3, ...
    'Parent',hp_SD, 'Tag','wave',...
    'Units','normalized', 'Position', [0.025 0.65 0.3 0.325],...
    'Callback',{@selectEventImage,mainFig}, 'Enable','on');

% button group to select type of event
h_bg_eventSelection = uibuttongroup('Visible','on', ...
    'FontUnits','normalized', 'FontSize',0.15,...
    'FontWeight','bold', 'Title','type of event:',...
    'Parent',hp_SD, 'Units','normalized', ...
    'Position', [0.025 0.025 0.3 0.6],...
    'SelectionChangedFcn','');
% Create radio buttons in the button group.
h_b1 = uicontrol(h_bg_eventSelection, 'Style','radiobutton', ...
    'Units','normalized',...
    'String','wave', 'FontUnits','normalized', 'FontSize',0.9,...
    'Position',[0.05 0.775 0.9 0.2],...
    'HandleVisibility','off');
h_b2 = uicontrol(h_bg_eventSelection, 'Style','radiobutton', ...
    'Units','normalized',...
    'String','transient (triggered)', ...
    'FontUnits','normalized', 'FontSize',0.9,...
    'Position',[0.05 0.525 0.9 0.2],...
    'HandleVisibility','off');
h_b3 = uicontrol(h_bg_eventSelection, 'Style','radiobutton', ...
    'Units','normalized',...
    'String','transient (spon.)', ...
    'FontUnits','normalized', 'FontSize',0.9,...
    'Position',[0.05 0.275 0.9 0.2],...
    'HandleVisibility','off');
h_b4 = uicontrol(h_bg_eventSelection, 'Style','radiobutton', ...
    'Units','normalized',...
    'String','caffeine', ...
    'FontUnits','normalized', 'FontSize',0.9,...
    'Position',[0.05 0.025 0.9 0.2],...
    'HandleVisibility','off');

% table of selected events
h_table_eventsImgs = uitable('Data',repmat({'----'},[5,1]), ...
    'ColumnEditable',true,...
    'ColumnName',{'<HTML> <font size="5"> <b> name of ROI </b> </HTML>'},...
    'ColumnWidth','auto', 'Parent',hp_SD,...
    'Units','normalized', ...
    'Position',[0.35 0.025 0.425 0.9],...
    'FontUnits','normalized', 'FontSize',0.10,...
    'CellEditCallback',{@tableEditEvents,mainFig});
h_table_eventsImgs.Units = 'pixels';
colWidth = h_table_eventsImgs.Position(3) - 50;
try
    h_table_eventsImgs.ColumnWidth = {colWidth};
catch
    h_table_eventsImgs.ColumnWidth = {h_table_eventsImgs.Position(3)};
end
% analyze events
h_pb_eventsAnalysis = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> analyze selected <br> events <html>',...
    'FontUnits','normalized', 'FontSize',0.225, ...
    'Parent',hp_SD, 'FontWeight','bold',...
    'Units','normalized', ...
    'Position', [0.8 0.25 0.175 0.5],...
    'Enable','on', 'Tag','calcParamsOfEvents', ...
    'Callback',{@selectedEventsAnalysisWindow,mainFig});

hObjs.h_push_eventImgROI = h_push_eventImgROI;
hObjs.h_bg_eventSelection = h_bg_eventSelection;
hObjs.h_b1_eventSelection = h_b1;
hObjs.h_b2_eventSelection = h_b2;
hObjs.h_b3_eventSelection = h_b3;
hObjs.h_b4_eventSelection = h_b4;
hObjs.h_table_eventsImgs = h_table_eventsImgs;
hObjs.h_pb_eventsAnalysis = h_pb_eventsAnalysis;


%% preview and save pushbutton 
hp_Analysis = uipanel('Title','save', 'Parent',mainFig,...
    'Units','normalized', ...
    'Position',[0.75 0.09 0.1 0.21],...
    'FontUnits','normalized', ...
    'FontSize',0.1,'FontWeight','bold');

h_pb_saveAnalysis = uicontrol('Style','pushbutton',...
    'String','save analysis',...
    'FontUnits','normalized', 'FontSize',0.45, ...
    'FontWeight','bold', 'Units','normalized',...
    'Parent',hp_Analysis, 'Position', [0.1 0.7 0.8 0.25],...
    'Callback',{@saveAnalysis,mainFig});

check_saveEventsFigs = uicontrol('Style','checkbox', ...
    'Parent',hp_Analysis,...
    'FontUnits','normalized', 'Value',1, ...
    'Units','normalized', 'FontSize',0.8,...
    'String','save individual events', ...
    'Position',[0.025 0.55 0.95 0.1]);

check_saveProfsAndFits = uicontrol('Style','checkbox', ...
    'Parent',hp_Analysis,...
    'FontUnits','normalized', 'Value',1, ...
    'Units','normalized', 'FontSize',0.8,...
    'String','save profs & fits', 'Position',[0.025 0.4 0.95 0.1]);

check_addTemperature = uicontrol('Style','checkbox', ...
    'Parent',hp_Analysis,...
    'FontUnits','normalized', 'Value',0, ...
    'Units','normalized', 'FontSize',0.8,...
    'String','add temperature record', ...
    'Position',[0.025 0.275 0.95 0.1]);

uicontrol('Style','text', ...
    'FontUnits','normalized',...
    'Parent',hp_Analysis, 'Units','normalized', 'FontSize',0.4,...
    'Position',[0.025 0.05 0.6 0.2], ...
    'String',{'output .pdf file';'resolution'});
h_edit_res = uicontrol('Style', 'edit', ...
    'String','300', 'Parent',hp_Analysis, ...
    'Units','normalized', 'FontUnits','normalized', ...
    'Position',[0.625 0.075 0.35 0.15], ...
    'Units','normalized', 'FontSize',0.6, ...
    'TooltipString','resolution of output pdf file',...
    'Callback','');

hObjs.h_pb_saveAnalysis = h_pb_saveAnalysis;
hObjs.check_saveProfsAndFits = check_saveProfsAndFits;
hObjs.check_addTemperature = check_addTemperature;
hObjs.check_saveEventsFigs = check_saveEventsFigs;
hObjs.h_edit_res = h_edit_res;


%% re-analyze
h_pb_reAnalyze = uicontrol('Style','pushbutton',...
    'String',['<html> <p align="center"> re-analyze <br/>', ...
              '<i> <span style="font-size:50%">(load image & .xls files)</span></i> <html>'],...
    'FontUnits','normalized','FontSize',0.3,'FontWeight','bold',...
    'Parent',mainFig,'Units','normalized','Position', [0.86 0.215 0.075 0.075],...
    'Callback',{@reAnalyze,mainFig},'Enable','on');

hObjs.h_pb_reAnalyze = h_pb_reAnalyze;


%% load temperature recording
hp_tempRec = uipanel('Title','temperature rec.', 'Parent',mainFig,...
    'Units','normalized', ...
    'Position',[0.86 0.09 0.075 0.12],...
    'FontUnits','normalized', ...
    'FontSize',0.13,'FontWeight','bold');

h_pb_selectTempDir = uicontrol('Style','pushbutton',...
    'String',['<html> <p align="center"> select dir <br/>', ...
              '<i> <span style="font-size:50%">(both temperature and <br/>', ...
              'stimulus record)</span></i> <html>'], ...
    'FontUnits','normalized', 'FontSize',0.3, ...
    'FontWeight','bold', 'Units','normalized',...
    'Parent',hp_tempRec, 'Position', [0.05 0.35 0.9 0.6],...
    'Callback',{@getTemperatureRecord,mainFig});

uicontrol('Style','text', ...
    'FontUnits','normalized',...
    'Parent',hp_tempRec, 'Units','normalized', 'FontSize',0.4,...
    'Position',[0.025 0.025 0.6 0.35], ...
    'HorizontalAlignment','center', ...
    'String',{'trigger', 'delay (s):'});
h_edit_triggerDelay = uicontrol('Style', 'edit', ...
    'String','26.5', 'Parent',hp_tempRec, ...
    'Units','normalized', 'FontUnits','normalized', ...
    'Position',[0.625 0.025 0.35 0.3], ...
    'Units','normalized', 'FontSize',0.55, ...
    'TooltipString','delay used in protocol to trigger microscope recording in (s)',...
    'Callback','');

hObjs.h_pb_selectTempDir = h_pb_selectTempDir;
hObjs.h_edit_triggerDelay = h_edit_triggerDelay;

%%

% save handles of objects
setappdata(mainFig,'hObjs',hObjs);


