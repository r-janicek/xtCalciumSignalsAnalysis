function createMaskWithROI(hO,~,analysisFig)

% data
hObjsA = getappdata(analysisFig,'hObjsA');
selectedEvent = getappdata(analysisFig,'selectedEvent');
t = selectedEvent.ROIdata.dataROIs.t;
x = selectedEvent.ROIdata.dataROIs.x;

% este pridat moznost volby miesta masky, mozno aj pre transient nejaky
% clonning, mozno dat iba nejaky noise okolo nuly
% switch selectedEvent.type
%     case 'wave'
%         currentAnalysis = selectedEvent.analysis.waveAnalysis;
%     case 'transient'
%         currentAnalysis = selectedEvent.analysis.transientAnalysis;
%     case 'caffeine'
%         currentAnalysis = selectedEvent.analysis.caffeineAnalysis;
%         
% end

% delete previous if any
delete(findall(findall(analysisFig,'Tag','tempROIaxes'),'Tag','maskROI'))
delete(findall(analysisFig,'Tag','tempROIaxes'))

% show whole image for better selection of mask
detImg = findall(hObjsA.ax_detEvent,'Type','image');
detImg.CData = selectedEvent.ROIdata.dataROIs.imgFN;

% create dummy axes only for selecting mask
ax_maskROI = axes('Parent',analysisFig,'Position',hObjsA.ax_detEvent.Position);
set(ax_maskROI,'XTick',[],'YTick',[],'Color','none','Tag','tempROIaxes')
% disable zooming, moving by mouse in axes
disableDefaultInteractivity(ax_maskROI)
% disableDefaultInteractivity(hObjsA.ax_detEvent) % maybe 

% create temporal transparent mask
tempMask = selectedEvent.detectedEvent.mask;
h_tempMaskImg = image(ax_maskROI,'CData',tempMask,...
    'CDataMapping','scaled','AlphaData',0.25);
set(ax_maskROI,'YLim',[0 size(tempMask,1)],'XLim',[0 size(tempMask,2)])
ax_maskROI.YAxis.Direction = 'reverse';

% save data
selectedEvent.detectedEvent.ax_maskROI = ax_maskROI;
selectedEvent.detectedEvent.tempMask = tempMask;
selectedEvent.detectedEvent.h_tempMaskImg = h_tempMaskImg;
setappdata(analysisFig,'selectedEvent',selectedEvent)

% draw ROI
% check mask proportions
hROI = drawrectangle(ax_maskROI, 'Color',[0 1 0], ...
    'Position',[1,1,3,3], ...
    'InteractionsAllowed','translate', ...
    'Tag', 'addMask');
% create context menu
cm = uicontextmenu(analysisFig);
m1 = uimenu(cm,'Text','add mask', ...
    'MenuSelectedFcn',{@setEventMaskCorrectionMode,hROI});
m2 = uimenu(cm,'Text','remove mask', ...
    'MenuSelectedFcn',{@setEventMaskCorrectionMode,hROI});
% assign it to mask correction ROI
hROI.ContextMenu = cm;

% notify during interaction with ROI
notify(hROI,'MovingROI') % 'MovingROI' ROIMoved
notify(hROI,'ROIClicked')

% listeners/callbacks
lh_moving = addlistener(hROI,'MovingROI',@updateMaskOfEventWithROI);
%lh_click = addlistener(hROI,'ROIClicked',@changeColorOfROI);

% set up functions to control size and position of ROI
analysisFig.WindowButtonMotionFcn = {@positionAndSizeOfcellMaskROI,hROI};
analysisFig.WindowScrollWheelFcn = {@positionAndSizeOfcellMaskROI,hROI};
if isempty(analysisFig.WindowScrollWheelFcn)
    % throw an error dlg
    opts = struct('WindowStyle', 'modal', 'Style', 'tex');
    errordlg({'Please turn off all interactivity in all axes!','e.g. zooming'}, ...
        'Error')
    keyboard
end

   

end

