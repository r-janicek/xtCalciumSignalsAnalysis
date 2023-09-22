function profilesFitWindow(~,~,mainFig)
% window to fit selected profiles with repetitive sparks
imgData = getappdata(mainFig,'imgData');
profileAnalysis = getappdata(mainFig,'profileAnalysis');

if ~isfield(profileAnalysis,'selectedROIs') || ...    
    isempty(profileAnalysis.selectedROIs)
    return
end

nProfiles = height(profileAnalysis.selectedROIs);

%% create figure for fitting  
fitFig = figure('Name','profile fit','units','normalized','outerposition',[0 0.2 1 0.8]);
set(mainFig, 'PaperPositionMode','auto','PaperOrientation',...
              'landscape','PaperType','A4','Tag','profilesFit');

h_txt_fitNum = uicontrol('Style','text',...
    'Parent',fitFig,'Units','normalized','Position', [0.1 0.96 0.8 0.03],...
    'FontUnits','normalized','FontSize',0.7,'FontWeight','bold',...
    'HorizontalAlignment','center','String','profile #1');

% create axes
ax_fit = axes('Parent',fitFig);
set(ax_fit,'Position',[0.03 0.6 0.94 0.35])
% set(get(ax_img,'Xlabel'),'String','t (ms)')
set(ax_fit,'XTick',[],'YGrid','on','FontSize',14)
set(get(ax_fit,'Ylabel'),'String','profile (F)','FontWeight','bold')
set(ax_fit,'buttondownfcn',{@mouseSetMaskFcn})

% plot first profile
prof = getSelectedProfileData(1,imgData,profileAnalysis,getappdata(mainFig,'analysisType'));


%  keyboard
%%%%%
% d = 5
% prof.t = d(:,1);
% prof.y = d(:,2);
% prof.baselineM = true(size(prof.y));
% prof.eventsM = false(size(prof.y));
% prof.eventsPeaks = {100,100};
% % 
% % %%%%%


line(prof.t,prof.y,'Parent',ax_fit,'Color','k','LineStyle','-','LineWidth',1,...
    'Tag','profile');
set(ax_fit,'XLim',[prof.t(1) prof.t(end)],'YLim',[min(prof.y)*0.95 max(prof.y)*1.05])
% show selected peaks
line([prof.eventsPeaks{:,2}],ones(size([prof.eventsPeaks{:,2}])).*ax_fit.YLim(2),...
    'Parent',ax_fit,'Color','r','LineStyle','none','LineWidth',1,...
    'Marker','.','MarkerSize',profileAnalysis.peaksCircleSz,'Tag','selectedPeaks');

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
    
ax_res = axes('Parent',fitFig);
set(ax_res,'Position',[0.03 0.39 0.94 0.2])
% set(get(ax_img_sparks,'Xlabel'),'String','t (ms)')
set(ax_res,'FontSize',14)
set(get(ax_res,'Ylabel'),'String','residuals (F)','FontWeight','bold')
set(get(ax_res,'Xlabel'),'String','t (ms)','FontWeight','bold')

hlink = linkprop([ax_fit,ax_res],'XLim');

hObjsFit.h_txt_fitNum = h_txt_fitNum;
hObjsFit.ax_fit = ax_fit;
hObjsFit.ax_res = ax_res;
hObjsFit.hlink = hlink;


%% create slider when image is too big in time direction
h_bg_sld_fit = uibuttongroup('Visible','on',...
    'Parent',fitFig,'Units','normalized','Position', [0.1 0.2275 0.1 0.08],...
                  'SelectionChangedFcn',{@bg_sld_selection});
                     
% Create three radio buttons in the button group.
h_b1 = uicontrol(h_bg_sld_fit,'Style','radiobutton','Units','normalized',...
                  'String','whole profile','FontUnits','normalized','FontSize',0.75,...
                  'Position',[0.05 0.55 0.9 0.35],...
                  'HandleVisibility','off');
              
h_b2 = uicontrol(h_bg_sld_fit,'Style','radiobutton','Units','normalized',...
                  'String','#s part','FontUnits','normalized','FontSize',0.75,...
                  'Position',[0.05 0.1 0.9 0.35],...
                  'HandleVisibility','off');

h_edit_imgSld_w_fit = uicontrol('Style','edit','String','1',...
                    'FontUnits','normalized','FontSize',0.75,...
                    'Parent',fitFig,'Units','normalized',...
                    'Position',[0.1 0.19 0.02 0.03],...
                    'Callback', {@bg_sld_selection});
uicontrol('Style','text','Parent',fitFig,'Units','normalized',...
            'Position', [0.125 0.19 0.01 0.03],...
            'FontUnits','normalized','FontSize',0.75,'String','s');
                
h_img_sld_fit = uicontrol('Style', 'slider','Parent',fitFig,'FontUnits','normalized',...
                'Min',1,'Max',100,'Value',1,'Units','normalized',...
                'Position', [0.03 0.59 0.94 0.01],'BackgroundColor',[0.7 0.7 0.7],...
                'Callback', '','Enable','off'); % {@imgShow_slider}

hObjsFit.h_bg_sld_fit = h_bg_sld_fit;
hObjsFit.h_b1_sld_fit = h_b1;
hObjsFit.h_b2_sld_fit = h_b2;
hObjsFit.h_edit_imgSld_w_fit = h_edit_imgSld_w_fit;

hObjsFit.h_img_sld_fit = h_img_sld_fit;


%% baseline
hpBs = uipanel('Title','baseline','Parent',fitFig,'Units','normalized',...
    'Position',[0.25 0.025 0.15 0.3],...
    'Visible', 'on','FontUnits','normalized',...
    'FontSize',0.08,'FontWeight','bold');

uicontrol('Style','text',...
    'Parent',hpBs,'Units','normalized','Position', [0.01 0.875 0.54 0.1],...
    'FontUnits','normalized','FontSize',0.725,'FontWeight','normal',...
    'HorizontalAlignment','right','String','baseline fit function:');

popUpMenuBs = uicontrol('Style','popup',...
    'Parent',hpBs,'Units','normalized','Position', [0.55 0.875 0.45 0.1],...
    'FontUnits','normalized','FontSize',0.725,'Tag','popUpMenuBsFcn',...
    'String',{'spline','SmoothingSpline','exp1','exp2','stretchedExp','polynomial'},...
    'Callback',{@selectBsFitFcn,fitFig});

% add slider to control baseline percentile 
startVal = floor(sum(prof.baselineM)/numel(prof.baselineM)*100);
if startVal==100, startVal=99; end
if startVal==0, startVal=1; end
h_txt_sld_Bs = uicontrol('Style','text',...
    'Parent',hpBs,'Units','normalized','Position', [0.05 0.725 0.9 0.1],...
    'FontUnits','normalized','FontSize',0.6,'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'String',sprintf('mask of baseline as %d percentile',startVal));

% set up slider, percentiles 
sldMin = 1;
sldMax = 99; 

sldStep = 1/(sldMax-sldMin);
sld_Bs = uicontrol('Style', 'slider','Parent',hpBs,'FontUnits','normalized',...
    'Min',sldMin,'Max',sldMax,'Value',startVal,'SliderStep',[sldStep 10*sldStep],...
    'Units','normalized','Position', [0.05 0.625 0.9 0.1],...
    'Callback', {@bsMask_slider,fitFig},'Enable','on');

h_txt_paramFitBs1 = uicontrol('Style','text',...
    'Parent',hpBs,'Units','normalized','Position', [0.4 0.4865 0.375 0.1],...
    'FontUnits','normalized','FontSize',0.725,'FontWeight','normal',...
    'HorizontalAlignment','right','String','fit param #1:');

h_edit_paramFitBs1 = uicontrol('Style','edit','String','nan',...
    'FontUnits','normalized','Parent',hpBs,'Units','normalized',...
    'FontSize',0.5,'Position', [0.8 0.475 0.15 0.125],...
    'Callback', '');

h_txt_paramFitBs2 = uicontrol('Style','text',...
    'Parent',hpBs,'Units','normalized','Position', [0.4 0.3365 0.375 0.1],...
    'FontUnits','normalized','FontSize',0.725,'FontWeight','normal',...
    'HorizontalAlignment','right','String','fit param #2:');

h_edit_paramFitBs2 = uicontrol('Style','edit','String','nan',...
    'FontUnits','normalized','Parent',hpBs,'Units','normalized',...
    'FontSize',0.5,'Position', [0.8 0.325 0.15 0.125],...
    'Callback', '');

h_pb_fitBs = uicontrol('Style','pushbutton','String','<html> <p align="center"> fit of <br> baseline <html>',...
    'FontUnits','normalized','Parent',hpBs,'Units','normalized','FontSize',0.3,...
    'Position', [0.05 0.325 0.3 0.275],'Callback', {@fitBaseline,fitFig});

h_pb_normProf = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> normalize <br> (&#916F/F0) <html>','FontUnits','normalized',...
    'FontSize',0.3,'FontWeight','bold',...
    'Parent',hpBs,'Units','normalized','Position', [0.55 0.01 0.4 0.275],...
    'Callback',{@normProfileFit,fitFig},'Enable','on');

h_pb_undoNormProf = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> undo <br> normalization <html>','FontUnits','normalized',...
    'FontSize',0.3,'FontWeight','normal',...
    'Parent',hpBs,'Units','normalized','Position', [0.05 0.01 0.4 0.275],...
    'Callback',{@undoNormProfileFit,fitFig},'Enable','on');

hObjsFit.popUpMenuBs = popUpMenuBs;
hObjsFit.h_txt_paramFitBs1 = h_txt_paramFitBs1;
hObjsFit.h_edit_paramFitBs1 = h_edit_paramFitBs1;
hObjsFit.h_txt_paramFitBs2 = h_txt_paramFitBs2;
hObjsFit.h_edit_paramFitBs2 = h_edit_paramFitBs2;
hObjsFit.h_pb_fitBs = h_pb_fitBs;
hObjsFit.h_pb_normProf = h_pb_normProf;
hObjsFit.h_pb_undoNormProf = h_pb_undoNormProf;
hObjsFit.h_txt_sld_Bs = h_txt_sld_Bs;
hObjsFit.sld_Bs = sld_Bs;

% set beggining
switch getappdata(mainFig,'analysisType')
    case 'spark recovery ryanodine'
        set(hObjsFit.popUpMenuBs,'Value',6)
        set(hObjsFit.h_txt_paramFitBs1,'String','order <1,9>')
        set(hObjsFit.h_edit_paramFitBs1,'String',num2str(4),'Enable','on')
        set(hObjsFit.h_txt_paramFitBs2,'String','spline order')
        set(hObjsFit.h_edit_paramFitBs2,'String','nan','Enable','off')
    otherwise
        set(hObjsFit.popUpMenuBs,'Value',6)
        set(hObjsFit.h_txt_paramFitBs1,'String','order <1,9>')
        set(hObjsFit.h_edit_paramFitBs1,'String',num2str(8),'Enable','on')
        set(hObjsFit.h_txt_paramFitBs2,'String','spline order')
        set(hObjsFit.h_edit_paramFitBs2,'String','nan','Enable','off')

%         set(hObjsFit.h_txt_paramFitBs1,'String','number of knots') 
%         set(hObjsFit.h_edit_paramFitBs1,'String',num2str(ceil(imgData.t(end)/200)),'Enable','on')
%         set(hObjsFit.h_txt_paramFitBs2,'String','order <2,10>')
%         set(hObjsFit.h_edit_paramFitBs2,'String',num2str(3),'Enable','on')
        
        prof.baselineM(1) = 1;
        prof.baselineM(end) = 1;
  
end

%% mask of events
maskButtonGroup = uibuttongroup(fitFig,'Units','normalized',...
                  'Title','mask','FontUnits','normalized',...
                  'FontSize',0.2,'FontWeight','bold',...
                  'Position', [0.41 0.18 0.1 0.12]);
            
% Create three radio buttons in the button group.
r1 = uicontrol('Style','radiobutton','Parent',maskButtonGroup,...
        'FontUnits','normalized','Value',0,'Units','normalized','FontSize',0.5,...
        'String','set mask of baseline','Position',[0.05 0.5 0.9 0.4]);
              
% Create three radio buttons in the button group.
r2 = uicontrol('Style','radiobutton','Parent',maskButtonGroup,...
        'FontUnits','normalized','Value',0,'Units','normalized','FontSize',0.5,...
        'String','set mask of events','Position',[0.05 0.05 0.9 0.4],'Enable','off');

hObjsFit.maskButtonGroup = maskButtonGroup;
hObjsFit.rbutton1 = r1;
hObjsFit.rbutton2 = r2;


%% select event function panel
hpEventFcn = uipanel('Title','event fit function','Parent',fitFig,'Units','normalized',...
    'Position',[0.41 0.05 0.1 0.12],...
    'Visible', 'on','FontUnits','normalized',...
    'FontSize',0.2,'FontWeight','bold');

popUpMenuEventsFcn = uicontrol('Style','popup',...
    'Parent',hpEventFcn,'Units','normalized','Position', [0.01 0.6 0.98 0.25],...
    'FontUnits','normalized','FontSize',0.85,'Tag','popUpMenuBsFcn','Value',4,...
    'String',{'expModGauss', '1expR1expD', '1expR2expD', 'CaSpikeFun'});

check_signOfPeak = uicontrol('Style', 'checkbox','Parent',hpEventFcn,...
    'FontUnits','normalized','Value',1,'Units','normalized','FontSize',0.85,...
    'String','positive peaks','Position', [0.01 0.1 0.98 0.25]);          
          
hObjsFit.popUpMenuEventsFcn = popUpMenuEventsFcn;
hObjsFit.check_signOfPeak = check_signOfPeak;


%% fit
hpFit = uipanel('Title','fit profile','Parent',fitFig,'Units','normalized',...
    'Position',[0.52 0.05 0.15 0.25],...
    'Visible', 'on','FontUnits','normalized',...
    'FontSize',0.1,'FontWeight','bold');

h_txt_fitNumSld = uicontrol('Style','text',...
    'Parent',hpFit,'Units','normalized','Position', [0.05 0.85 0.9 0.15],...
    'FontUnits','normalized','FontSize',0.5,'FontWeight','bold',...
    'HorizontalAlignment','center','String','profile #1');

% set up slider 
sldMin = 1;
sldMax = nProfiles; % nProfiles

sldStep = 1/(sldMax-sldMin);
if sldStep>=1, sldStep=1; end

sld_fitNum = uicontrol('Style', 'slider','Parent',hpFit,'FontUnits','normalized',...
    'Min',1,'Max',nProfiles,'Value',1,'SliderStep',[sldStep sldStep],...
    'Units','normalized','Position', [0.05 0.7 0.9 0.15],...
    'Callback', {@fitNum_slider,fitFig},'Enable','on');

if sldStep == 1 && sldMin==sldMax
    set(sld_fitNum,'Enable','off')
end
%
h_pb_fit = uicontrol('Style','pushbutton','String','fit profile',...
    'FontUnits','normalized','Parent',hpFit,'Units','normalized','FontSize',0.5,...
    'Position', [0.3 0.475 0.4 0.2],'Callback', {@fitProfile,fitFig});

check_fitWeights = uicontrol('Style', 'checkbox','Parent',hpFit,...
    'FontUnits','normalized','Value',0,'Units','normalized','FontSize',0.85,...
    'String','weights','Position', [0.725 0.575 0.25 0.1]);        

h_pb_acceptFit = uicontrol('Style', 'pushbutton',...
    'String','accept fit','FontUnits','normalized',...
    'FontSize',0.5,'FontWeight','bold',...
    'Parent',hpFit,'Units','normalized','Position', [0.3 0.25 0.4 0.2],...
    'Callback',{@acceptFit,fitFig},'Enable','on');

h_pb_repairFit = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> repair of<br> accepted fit <html>','FontUnits','normalized',...
    'FontSize',0.3,'FontWeight','normal',...
    'Parent',hpFit,'Units','normalized','Position', [0.3 0.025 0.4 0.2],...
    'Callback',{@repairAcceptedFit,fitFig},'Enable','on');

hObjsFit.h_txt_fitNumSld = h_txt_fitNumSld;
hObjsFit.sld_fitNum = sld_fitNum;
hObjsFit.h_pb_fit = h_pb_fit;
hObjsFit.check_fitWeights = check_fitWeights;
hObjsFit.h_pb_acceptFit = h_pb_acceptFit;
hObjsFit.h_pb_repairFit = h_pb_repairFit;

% 
h_pb_stopFiting = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> stop <br> current fitting <html>','FontUnits','normalized',...
    'FontSize',0.25,'FontWeight','normal',...
    'Parent',fitFig,'Units','normalized','Position', [0.68 0.2 0.07 0.08],...
    'Callback',{@stopFitting,fitFig},'Enable','on');

%
h_pb_closeFitWin = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> close <br> fit window <html>','FontUnits','normalized',...
    'FontSize',0.25,'FontWeight','normal',...
    'Parent',fitFig,'Units','normalized','Position', [0.68 0.075 0.07 0.1],...
    'Callback',{@closeFitWin,fitFig},'Enable','on');

hObjsFit.h_pb_stopFiting = h_pb_stopFiting;
setappdata(fitFig,'stopFiting',false)

%% save data
setappdata(fitFig,'hObjsFit',hObjsFit)
setappdata(fitFig,'imgData',imgData)
setappdata(fitFig,'profileAnalysis',profileAnalysis)
setappdata(fitFig,'mainFig',mainFig)


% save selected profile data
setappdata(fitFig,'selectedProf',prof) 

% do fitting of baseline
fitBaseline([],[],fitFig)

end

