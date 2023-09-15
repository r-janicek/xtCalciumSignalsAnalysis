function selectedEventsAnalysisWindowOLD(~,~,mainFig)

% get data
imgData = getappdata(mainFig,'imgData');

if ~isfield(getappdata(mainFig),'selectedROIs') || isempty(getappdata(mainFig,'selectedROIs'))
    return
end

selectedROIs = getappdata(mainFig,'selectedROIs');

nROIs = height(selectedROIs);

%% create figure for analysis  
analysisFig = figure('Name','analysis of selected ROIs','units','normalized','outerposition',[0 0.2 1 0.8]);
set(mainFig, 'PaperPositionMode','auto','PaperOrientation',...
              'landscape','PaperType','A4','Tag','analysis of selected ROIs');

h_txt_fitNum = uicontrol('Style','text',...
    'Parent',analysisFig,'Units','normalized','Position', [0.1 0.96 0.8 0.03],...
    'FontUnits','normalized','FontSize',0.7,'FontWeight','bold',...
    'HorizontalAlignment','center','String',sprintf('analysis of: %s',selectedROIs.roiName{1}));

% create axes
ax_orgImg = axes('Parent',analysisFig);
set(ax_orgImg,'Position',[0.04 0.55 0.28 0.35])
% show first event
image(selectedROIs.dataROIs(1).imgFN,'YData',[min(selectedROIs.dataROIs(1).x) max(selectedROIs.dataROIs(1).x)],...
    'XData',[min(selectedROIs.dataROIs(1).t) max(selectedROIs.dataROIs(1).t)],...
    'CDataMapping','scaled','Parent',ax_orgImg);
set(ax_orgImg,'XTick',[],'YGrid','off','FontSize',14)
set(get(ax_orgImg,'Ylabel'),'String','x (\mum)','FontWeight','bold')
set(get(ax_orgImg,'title'),'String','image filtered and normalized','FontWeight','bold')

ax_orgImgProf = axes('Parent',analysisFig);
set(ax_orgImgProf,'Position',[0.04 0.34 0.28 0.2])
set(ax_orgImgProf,'YGrid','on','FontSize',14)
set(get(ax_orgImgProf,'Ylabel'),'String','fluorescence (\DeltaF/F0)','FontWeight','bold')
set(get(ax_orgImgProf,'Xlabel'),'String','t (ms)','FontWeight','bold')
% profile of 1. event 
line(selectedROIs.dataROIs(1).t,mean(selectedROIs.dataROIs(1).imgFN,1),...
    'Parent',ax_orgImgProf,'Color','k','LineStyle','-','LineWidth',1,...
    'Tag','profile');
set(ax_orgImgProf,'XLim',[min(selectedROIs.dataROIs(1).t) max(selectedROIs.dataROIs(1).t)],...
    'YLim',[min(mean(selectedROIs.dataROIs(1).imgFN,1))*0.95 max(mean(selectedROIs.dataROIs(1).imgFN,1))*1.05])  

ax_orgImgRes = axes('Parent',analysisFig);
set(ax_orgImgRes,'Position',[0.04 0.13 0.28 0.15])
set(ax_orgImgRes,'YGrid','on','FontSize',14)
set(get(ax_orgImgRes,'Ylabel'),'String','residuals','FontWeight','bold')
set(get(ax_orgImgRes,'Xlabel'),'String','t (ms)','FontWeight','bold')


ax_detEvent = axes('Parent',analysisFig);
set(ax_detEvent,'Position',[0.36 0.55 0.28 0.35])
% show dummy image of size of first image
image(zeros(size(selectedROIs.dataROIs(1).imgFN)),'YData',[min(selectedROIs.dataROIs(1).x) max(selectedROIs.dataROIs(1).x)],...
    'XData',[min(selectedROIs.dataROIs(1).t) max(selectedROIs.dataROIs(1).t)],...
    'CDataMapping','scaled','Parent',ax_detEvent);
set(ax_detEvent,'YGrid','off','FontSize',14)
set(get(ax_detEvent,'Ylabel'),'String','x (\mum)','FontWeight','bold')
set(get(ax_detEvent,'Xlabel'),'String','t (ms)','FontWeight','bold')
set(get(ax_detEvent,'title'),'String','detected event','FontWeight','bold')

ax_deskewed = axes('Parent',analysisFig);
set(ax_deskewed,'Position',[0.68 0.55 0.28 0.35])
% show dummy image of size of first image
image(zeros(size(selectedROIs.dataROIs(1).imgFN)),'YData',[min(selectedROIs.dataROIs(1).x) max(selectedROIs.dataROIs(1).x)],...
    'XData',[min(selectedROIs.dataROIs(1).t) max(selectedROIs.dataROIs(1).t)],...
    'CDataMapping','scaled','Parent',ax_deskewed);
set(ax_deskewed,'XTick',[],'YGrid','off','FontSize',14)
set(get(ax_deskewed,'Ylabel'),'String','x (\mum)','FontWeight','bold')
set(get(ax_deskewed,'title'),'String','deskewed event','FontWeight','bold')

ax_deskewedProf = axes('Parent',analysisFig);
set(ax_deskewedProf,'Position',[0.68 0.34 0.28 0.2])
set(ax_deskewedProf,'FontSize',14,'YGrid','on')
set(get(ax_deskewedProf,'Ylabel'),'String','fluorescence (\DeltaF/F0)','FontWeight','bold')
set(get(ax_deskewedProf,'Xlabel'),'String','t (ms)','FontWeight','bold')

ax_deskewedRes = axes('Parent',analysisFig);
set(ax_deskewedRes,'Position',[0.68 0.13 0.28 0.15])
set(ax_deskewedRes,'YGrid','on','FontSize',14)
set(get(ax_deskewedRes,'Ylabel'),'String','residuals','FontWeight','bold')
set(get(ax_deskewedRes,'Xlabel'),'String','t (ms)','FontWeight','bold')

hlink_org = linkprop([ax_orgImg,ax_orgImgProf],'XLim');
hlink_des = linkprop([ax_deskewed,ax_deskewedProf],'XLim');

hObjsA.h_txt_fitNum = h_txt_fitNum;
hObjsA.ax_orgImg = ax_orgImg;
hObjsA.ax_orgImgProf = ax_orgImgProf;
hObjsA.ax_orgImgRes = ax_orgImgRes;
hObjsA.ax_detEvent = ax_detEvent;
hObjsA.ax_deskewed = ax_deskewed;
hObjsA.ax_deskewedProf = ax_deskewedProf;
hObjsA.ax_deskewedRes = ax_deskewedRes;
hObjsA.hlink_org = hlink_org;
hObjsA.hlink_des = hlink_des;

% show photolytic pulses positions if there are any 
if isfield(imgData,'s_TPP')
    
    s_TPP  = imgData.s_TPP;
    e_TPP = imgData.e_TPP;
    
    if ~isempty(s_TPP)
        
        Xpp = [prof.t(s_TPP); 
               prof.t(e_TPP);
               prof.t(e_TPP);
               prof.t(s_TPP)];    
           
        Ypp = [ones(size(s_TPP)).*min(ax_fit.YLim);
               ones(size(s_TPP)).*min(ax_fit.YLim);
               ones(size(s_TPP)).*max(ax_fit.YLim);
               ones(size(s_TPP)).*max(ax_fit.YLim)];
        
        patch('XData',Xpp, 'YData',Ypp, 'Parent',ax_fit, 'Tag','photolyticPulses',...
            'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
    end
    
end


%% select event function panel
hpEventFcn = uipanel('Title','event fit function','Parent',analysisFig,'Units','normalized',...
    'Position',[0.36 0.22 0.12 0.25],...
    'Visible', 'on','FontUnits','normalized',...
    'FontSize',0.1,'FontWeight','bold');

popUpMenuEventsFcn = uicontrol('Style','popup',...
    'Parent',hpEventFcn,'Units','normalized','Position', [0.05 0.85 0.90 0.1],...
    'FontUnits','normalized','FontSize',0.85,'Tag','popUpMenuEvtFcn',...
    'String',{'spline', '1expR1expD', 'CaSpikeFun'},...
    'Callback',{@selectFitFcn,analysisFig});

h_txt_paramFit1 = uicontrol('Style','text',...
    'Parent',hpEventFcn,'Units','normalized','Position', [0.01 0.6 0.3 0.2],...
    'FontUnits','normalized','FontSize',0.4,'FontWeight','normal',...
    'HorizontalAlignment','center','String',{'fit param'; '#1:'});

h_edit_paramFit1 = uicontrol('Style','edit','String','nan',...
    'FontUnits','normalized','Parent',hpEventFcn,'Units','normalized',...
    'FontSize',0.6,'Position', [0.32 0.6 0.17 0.175],...
    'Callback', '');

h_txt_paramFit2 = uicontrol('Style','text',...
    'Parent',hpEventFcn,'Units','normalized','Position', [0.5 0.6 0.3 0.2],...
    'FontUnits','normalized','FontSize',0.4,'FontWeight','normal',...
    'HorizontalAlignment','center','String',{'fit param'; '#2:'});

h_edit_paramFit2 = uicontrol('Style','edit','String','nan',...
    'FontUnits','normalized','Parent',hpEventFcn,'Units','normalized',...
    'FontSize',0.6,'Position', [0.82 0.6 0.17 0.175],...
    'Callback', '');

check_pieceWise = uicontrol('Style', 'checkbox','Parent',hpEventFcn,...
    'FontUnits','normalized','Value',1,'Units','normalized','FontSize',0.85,...
    'String','piecewise function','Position', [0.1 0.45 0.9 0.1],...
    'Callback', {@selectFitFcn,analysisFig});  

uicontrol('Style','text',...
    'Parent',hpEventFcn,'Units','normalized','Position', [0.05 0.275 0.6 0.1],...
    'FontUnits','normalized','FontSize',0.85,'FontWeight','normal',...
    'HorizontalAlignment','center','String','number of peaks:');

h_edit_Npeaks = uicontrol('Style','edit','String','1',...
    'FontUnits','normalized','Parent',hpEventFcn,'Units','normalized',...
    'FontSize',0.6,'Position', [0.65 0.225 0.15 0.175],...
    'Callback', {@selectFitFcn,analysisFig});

check_signOfPeak = uicontrol('Style', 'checkbox','Parent',hpEventFcn,...
    'FontUnits','normalized','Value',1,'Units','normalized','FontSize',0.85,...
    'String','positive peaks','Position', [0.1 0.075 0.9 0.1]);    
     
hObjsA.popUpMenuEventsFcn = popUpMenuEventsFcn;
hObjsA.h_edit_Npeaks = h_edit_Npeaks;
hObjsA.check_signOfPeak = check_signOfPeak;
hObjsA.check_pieceWise = check_pieceWise;
hObjsA.h_txt_paramFit1 = h_txt_paramFit1;
hObjsA.h_edit_paramFit1 = h_edit_paramFit1;
hObjsA.h_txt_paramFit2 = h_txt_paramFit2;
hObjsA.h_edit_paramFit2 = h_edit_paramFit2;

% set beggining
setappdata(analysisFig,'hObjsA',hObjsA)
selectFitFcn([],[],analysisFig)


%% select baseline function panel
hpBaselineFcn = uipanel('Title','baseline fit function','Parent',analysisFig,'Units','normalized',...
    'Position',[0.36 0.14 0.12 0.07],...
    'Visible', 'on','FontUnits','normalized',...
    'FontSize',0.3,'FontWeight','bold');

popUpMenuBaselineFcn = uicontrol('Style','popup',...
    'Parent',hpBaselineFcn,'Units','normalized','Position', [0.05 0.05 0.90 0.9],...
    'FontUnits','normalized','FontSize',0.55,'Tag','popUpMenuBsFcn',...
    'String',{'exp1', 'exp2', 'stretchedExp', 'poly1', 'poly2', 'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'poly9'},...
    'Callback','','Value',8);

hObjsA.popUpMenuBaselineFcn = popUpMenuBaselineFcn;


%% smooth profile
hpSmooth = uipanel('Title','smooth profile','Parent',analysisFig,'Units','normalized',...
    'Position',[0.36 0.04 0.12 0.09],...
    'Visible', 'on','FontUnits','normalized',...
    'FontSize',0.27778,'FontWeight','bold');

h_smooth_txt = uicontrol('Style','text',...
    'Parent',hpSmooth,'Units','normalized','FontUnits','normalized','FontSize',0.33,...
    'Position', [0.01 0.05 0.64 0.9],'String',{'span of smoothing:';'(ms)\loess'});

h_smooth_edit = uicontrol('Style','edit','String','100',...
    'FontUnits','normalized','FontSize',0.6,...
    'Parent',hpSmooth,'Units','normalized','Position', [0.7 0.15 0.25 0.7],...
    'Callback',{@smoothProf,analysisFig});

hObjsA.h_smooth_edit = h_smooth_edit;


%% fit
hpFit = uipanel('Title','fit event','Parent',analysisFig,'Units','normalized',...
    'Position',[0.49 0.22 0.15 0.25],...
    'Visible', 'on','FontUnits','normalized',...
    'FontSize',0.1,'FontWeight','bold');

h_txt_fitNumSld = uicontrol('Style','text',...
    'Parent',hpFit,'Units','normalized','Position', [0.05 0.85 0.9 0.15],...
    'FontUnits','normalized','FontSize',0.5,'FontWeight','bold',...
    'HorizontalAlignment','center','String','event #1');

% set up slider 
sldMin = 1;
sldMax = nROIs; % nProfiles

sldStep = 1/(sldMax-sldMin);
if sldStep>=1, sldStep=1; end

sld_fitNum = uicontrol('Style', 'slider','Parent',hpFit,'FontUnits','normalized',...
    'Min',1,'Max',nROIs,'Value',1,'SliderStep',[sldStep sldStep],...
    'Units','normalized','Position', [0.05 0.7 0.9 0.15],...
    'Callback', {@fitNum_slider,analysisFig},'Enable','on');

if sldStep == 1 && sldMin==sldMax
    set(sld_fitNum,'Enable','off')
end

%
check_fitWeights = uicontrol('Style', 'checkbox','Parent',hpFit,...
    'FontUnits','normalized','Value',0,'Units','normalized','FontSize',0.85,...
    'String','weights','Position', [0.05 0.475 0.35 0.1]);        

h_pb_stopFiting = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> stop <br> current fitting <html>','FontUnits','normalized',...
    'FontSize',0.25,'FontWeight','normal',...
    'Parent',hpFit,'Units','normalized','Position', [0.05 0.05 0.35 0.3],...
    'Callback',{@stopFitting,analysisFig},'Enable','on');

h_pb_fit = uicontrol('Style','pushbutton','String','fit event',...
    'FontUnits','normalized','Parent',hpFit,'Units','normalized','FontSize',0.4,...
    'Position', [0.45 0.375 0.5 0.3],'Callback', {@fitSelectedEvent,analysisFig});

h_pb_acceptFit = uicontrol('Style', 'pushbutton',...
    'String','accept fit','FontUnits','normalized',...
    'FontSize',0.4,'FontWeight','bold',...
    'Parent',hpFit,'Units','normalized','Position', [0.45 0.05 0.5 0.3],...
    'Callback',{@acceptFit,analysisFig},'Enable','on');

%
h_pb_closeFitWin = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> close <br> fit window <html>','FontUnits','normalized',...
    'FontSize',0.25,'FontWeight','normal',...
    'Parent',analysisFig,'Units','normalized','Position', [0.5275 0.1 0.075 0.1],...
    'Callback',{@closeFitWin,analysisFig},'Enable','on');

hObjsA.h_txt_fitNumSld = h_txt_fitNumSld;
hObjsA.sld_fitNum = sld_fitNum;
hObjsA.h_pb_fit = h_pb_fit;
hObjsA.check_fitWeights = check_fitWeights;
hObjsA.h_pb_acceptFit = h_pb_acceptFit;
hObjsA.h_pb_stopFiting = h_pb_stopFiting;
setappdata(analysisFig,'stopFiting',false)


%% save data
setappdata(analysisFig,'hObjsA',hObjsA)
setappdata(analysisFig,'imgData',imgData)

% detect and analyze first event
outAnalysis = analyzeEvent(1,selectedROIs,imgData.pxSzT,imgData.pxSzX,analysisFig);

setappdata(analysisFig,'eventAnalysis',outAnalysis)
setappdata(analysisFig,'mainFig',mainFig)


end

