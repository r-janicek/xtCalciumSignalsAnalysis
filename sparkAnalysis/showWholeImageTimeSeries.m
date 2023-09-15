function showWholeImageTimeSeries(~,~,setUpFig,mainFig)

keyboard

% get data 
hObjs = getappdata(mainFig,'hObjs');
imgDataS = getappdata(setUpFig,'imgData');
tblData = getappdata(setUpFig,'h_tbl');

wSeriesFluoImgF = [];
wSeriesFluoImgR = [];
wSeriesTransImg = [];
t_wSeriesImgReal = [];
for i=1:numel(imgDataS)
    
    wSeriesFluoImgF = [wSeriesFluoImgF,imgDataS(i).imgDataXTfluoF]; % filtered
    wSeriesFluoImgR = [wSeriesFluoImgR,imgDataS(i).imgDataXTfluoR]; % raw
    wSeriesTransImg = [wSeriesTransImg,imgDataS(i).imgDataXTtrans]; % trans image
    t_img = imgDataS(i).t+tblData.Data{i,3}*1000;
    t_wSeriesImgReal = [t_wSeriesImgReal,t_img];
    
    szImg(i,:) = size(imgDataS(i).imgDataXTfluoF);    
end

t_px = (1:numel(t_wSeriesImgReal));

% create structure for whole series image, same as from openImg function, 
% and set up main figure window

% assign channels and set up popup menus
chNames = arrayfun(@(x) sprintf('ch #%d',x) ,(1:1), 'UniformOutput',0);
chNames = [chNames,{'none'}];

set(hObjs.popUpMenuChFluo,'String',chNames,'Value',1)
set(hObjs.popUpMenuChTrans,'String',chNames,'Value',2)

% plot image
ax_img = hObjs.ax_img;
cla(ax_img)

image(wSeriesFluoImgF,'YData',[1 size(wSeriesFluoImgF,1)],'XData',t_px,...
      'CDataMapping','scaled','Parent',ax_img);
set(ax_img,'XTick',[],'FontSize',14)
set(get(ax_img,'Ylabel'),'String','x (px)','FontWeight','bold')

% time profile 
prof_t = mean(wSeriesFluoImgF,1);

% clear axes
ax_img_sparks = hObjs.ax_img_sparks;
cla(ax_img_sparks)
% ax_img_sparks_2 = hObjs.ax_img_sparks_2;
% cla(ax_img_sparks_2)

ax_prof = hObjs.ax_prof;
plot(t_px,prof_t,'Parent',ax_prof,'Color','k');
set(ax_prof,'FontSize',14)
set(get(ax_prof,'Xlabel'),'String','t (px)','FontWeight','bold')
set(get(ax_prof,'Ylabel'),'String','Fluorescence (F)','FontWeight','bold')
xlim(ax_prof,[min(t_px) max(t_px)])

try
    ylim(ax_prof,[min(prof_t)*0.9 max(prof_t)])
catch
    ylim(ax_prof,[min(prof_t)*0.9 max(prof_t)+1])
end

% show patch object with drug application
m_ctrl = [tblData.Data{:,2}];

drugName = tblData.Data{find(m_ctrl==0,1,'first'),4};
drugConc = tblData.Data{find(m_ctrl==0,1,'first'),5};

yh_patch = abs(diff(ax_prof.YLim))*0.10;

szXimgCtrl = szImg(m_ctrl,2);
szXimgDrug = szImg(~m_ctrl,2);

cInc1 = 0.5/numel(szXimgCtrl);
cInc2 = 0.6/numel(szXimgDrug);

t_imgs = [tblData.Data{:,3}];
t_imgsCtrl = t_imgs(m_ctrl);
t_imgsDrug = t_imgs(~m_ctrl);

for i=1:numel(szXimgCtrl)
    
    ts = szXimgCtrl(i)*(i-1)+1;
    te = szXimgCtrl(i)*i;
 
    h_ctrlPatch = patch('XData',[ts te te ts],...
        'YData',[ax_prof.YLim(1) ax_prof.YLim(1) ax_prof.YLim(1)+yh_patch ax_prof.YLim(1)+yh_patch],...
        'EdgeColor','none','FaceColor','b','FaceAlpha',cInc1*i,'Parent',ax_prof);
    h_ctrlText = text(ts,ax_prof.YLim(1),sprintf('CTRL (img #%d, ts=%d s)',i,t_imgsCtrl(i)),'Color','w','FontSize',16,...
        'Parent',ax_prof,'VerticalAlignment','bottom','FontWeight','bold');
end

for i=1:numel(szXimgDrug)
    
    ts = szXimgDrug(i)*(i-1)+1+sum(szXimgCtrl);
    te = szXimgDrug(i)*i+sum(szXimgCtrl);
    
    h_drugPatch = patch('XData',[ts te te ts],...
        'YData',[ax_prof.YLim(1) ax_prof.YLim(1) ax_prof.YLim(1)+yh_patch ax_prof.YLim(1)+yh_patch],...
        'EdgeColor','none','FaceColor','r','FaceAlpha',cInc2*i,'Parent',ax_prof);
    h_drugText = text(ts,ax_prof.YLim(1),...
        sprintf('%s %d uM (img #%d, ts=%d s)',drugName,drugConc,i,t_imgsDrug(i)),'Color','k','FontSize',16,...
        'Parent',ax_prof,'VerticalAlignment','bottom','FontWeight','bold');
    
end

% link axes
linkaxes([ax_img,ax_img_sparks,ax_prof],'x')

% load positions of photolytic pulses, from trans channel 
if strcmp(getappdata(mainFig,'analysisType'),'spark recovery photolysis')

   [s_TPP,e_TPP,durOfTPP,TPP_delays,...
    d_TP_laser,posOfTPPinScanLine] = loadPhotolysisPositions(wSeriesFluoImgF,...
                                                             wSeriesTransImg,...
                                                             t_wSeriesImgReal,...
                                                             imgDataS.pxSzT,...
                                                             imgDataS.pxSzX,...
                                                             imgDataS.pos_P,...
                                                             imgDataS.pos_l,...
                                                             hObjs.ax_img,...
                                                             1);    
else
   
   TPP_delays = [nan;nan];
   posOfTPPinScanLine = [];
   durOfTPP = [];
   s_TPP = [];
   e_TPP = [];
   d_TP_laser = 0;
end

% maybe use percentile 1 and 99
z_max_img = max(wSeriesFluoImgF(:));
z_min_img = min(wSeriesFluoImgF(:));

% remove old image data
if isfield(getappdata(mainFig),'imgData')    
    rmappdata(mainFig,'imgData');    
end

% construct file name for series
seriesName = ['SeriesOfImgs_',strjoin(regexp([imgDataS.fileName],'\d*','match'),'_')];

% save image data
imgData = struct('filePath',imgDataS(1).filePath,...
                 'fileName',seriesName,...
                 'imgDataXT',wSeriesFluoImgR,...
                 'imgDataXY',{imgDataS(1).imgDataXY},...
                 'wholeImgXT',wSeriesFluoImgR,...
                 'wholeImgFluoXT',wSeriesFluoImgF,...   % image data, filtered
                 'imgDataXTfluoR',wSeriesFluoImgR,...  % raw image data
                 'imgDataXTfluoRN',wSeriesFluoImgR,... % raw image data, normalized
                 'imgDataXTfluoF',wSeriesFluoImgF,...   % image data, filtered
                 'imgDataXTfluoFN',wSeriesFluoImgF,...   % image data, filtered and normalized
                 'imgDataXTtrans',wSeriesTransImg,... % image data transmitted light
                 'pxSzT',imgDataS(1).pxSzT,...
                 'pxSzX',imgDataS(1).pxSzX,...
                 't',t_px,...
                 't_real',t_wSeriesImgReal,...
                 'm_ctrl',m_ctrl,...    
                 'szImgs',szImg,...                     % in px
                 't_imgs',t_imgs,...         % in s
                 'drugName',drugName,...
                 'drugConc',drugConc,...
                 'imgFiltersUsed',imgDataS(1).imgFiltersUsed,...
                 'z_max_img',z_max_img,...
                 'z_min_img',z_min_img,...
                 's_TPP',s_TPP,...
                 'e_TPP',e_TPP,...
                 'durOfTPP',durOfTPP,...
                 'TPP_delays',TPP_delays,...
                 'd_TP_laser',d_TP_laser,...
                 'TPPpointPos',[],...
                 'scanLinePos',imgDataS(1).scanLinePos,...
                 'posOfTPPinScanLine',posOfTPPinScanLine);

setappdata(mainFig,'imgData',imgData)
setappdata(mainFig,'seriesOfImgs',true)

% set up main window
set(hObjs.txt_name_img,'String',sprintf('%s',[imgData.filePath,char(seriesName)]));
set(hObjs.txt_name_img,'ToolTipString',sprintf('%s',[imgData.filePath,char(seriesName)]));
set(hObjs.h_table_profs,'Data',zeros(5,2))
% show trans image in small window
image(wSeriesTransImg,'CDataMapping','scaled','Parent',hObjs.h_ax_transCh,...
    'ButtonDownFcn',{@axesTransChFunction});
set(hObjs.h_ax_transCh,'XTick',[],'YTick',[],'PickableParts','all')
set(hObjs.h_ax_transCh,'ButtonDownFcn',{@axesTransChFunction})

% set normalization panel
set(hObjs.h_pb_norm,'Enable','on')
set(hObjs.h_pb_norm,'String','normalize data')
set(hObjs.h_pb_norm,'Callback',{@norm_data,mainFig})
set(hObjs.h_pb_norm,'FontWeight','normal')

% set blank panel
set(hObjs.h_push_BlankROI,'String','set ROI for blank','FontWeight','normal')
set(hObjs.h_push_BlankROI,'Callback', {@setBlankROI,mainFig})
set(hObjs.h_push_BlankROI,'Enable','on')

set(hObjs.h_edit_Blank,'String','nan')
set(hObjs.h_edit_Blank,'Enable','on')

if isempty(imgData.imgDataXY)    
    set(hObjs.h_pb_GetBlank,'Enable','off')    
else
    set(hObjs.h_pb_GetBlank,'Enable','on')   
end

% set up pushbutton for cropping
set(hObjs.h_pb_crop,'String','<html> <p align="center"> set ROI <br> for cropping <html>')
set(hObjs.h_pb_crop,'FontWeight','normal')


% set up panel for spark detection
% set(hObjs.check_visRect,'Value',1)
set(hObjs.txt_SpF,'String','# sp*100um-1*s-1')

if strcmp(get(hObjs.h_pb_crop,'String'),'crop')  
    
    set(hObjs.h_pb_crop,'String','set ROI for cropping')
    set(hObjs.h_pb_crop,'FontWeight','normal')
    
end
% set FWHM for spark detection, in pixels, 2 um
set(hObjs.h_edit_FWHM,'String',num2str(ceil(2/imgData.pxSzX)))


% set up panel for spark recovery
% set width of ROI for spark recovery to 10 um width in pixels (odd number)
n_px = round(10/imgData.pxSzX);
if mod(n_px,2)==0
   n_px = n_px+1; 
end

set(hObjs.h_edit_ROI,'String',num2str(n_px))


% remove old analysis
if isfield(getappdata(mainFig),'lastPressedPusbutton')
    rmappdata(mainFig,'lastPressedPusbutton');   
end

if isfield(getappdata(mainFig),'profileAnalysis')
    rmappdata(mainFig,'profileAnalysis');   
end


% restore image slider
imgSld_w = str2double(hObjs.h_edit_imgSld_w.String)*1000; % change from s to ms
set(hObjs.h_img_sld,'Value',1,'Max',ceil(imgData.t(end)/imgSld_w),'Enable','off')
% set XLim of axes
val_l = imgData.t(1);
val_u = imgData.t(end);
h_ax = [hObjs.ax_img,hObjs.ax_img_sparks,hObjs.ax_prof];
arrayfun(@(x) set(x,'XLim',[val_l val_u]),h_ax)
hObjs.h_bg_sld.SelectedObject = hObjs.h_b1_sld;

% remove listeners
delete(getappdata(mainFig,'propListener'))
delete(getappdata(mainFig,'editScaleListener'))

% remove main figure functions
set(mainFig,'WindowKeyPressFcn',[]);
set(mainFig,'WindowButtonUpFcn',[]);
set(mainFig,'WindowButtonMotionFcn',[]);


% close(setUpFig)

end

