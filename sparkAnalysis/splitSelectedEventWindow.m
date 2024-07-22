function splitSelectedEventWindow(mainFig, eventToSplit)
              
% create main figure and axes   
wFig = figure('Name','split event','units',...
    'normalized','outerposition',[0.2 0.3 0.6 0.7]);
set(wFig, 'PaperPositionMode','auto','PaperOrientation',...
              'landscape','PaperType','A4','Tag','sparkAnalysis');

% create axes
ax_img = axes('Parent',wFig);
set(ax_img,'Position',[0.075 0.4 0.9 0.575])
% set(get(ax_img,'Xlabel'),'String','t (ms)')
set(ax_img,'XTick',[],'FontSize',14)
set(get(ax_img,'Ylabel'),'String','x (pixels)','FontWeight','bold')
hObjsW.ax_img = ax_img;

% panel for watershed transform
hp_SD = uipanel('Title','split events','Parent',wFig,...
    'Position',[0.15 0.025 0.7 0.275],'FontUnits','normalized',...
    'FontSize',0.1,'FontWeight','bold');

h_push_doWatershed = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> do watershed <br> transform <html>',...
    'FontUnits','normalized','FontSize',0.25, ...
    'FontWeight','bold','Parent',hp_SD,...
    'Units','normalized','Position', [0.115 0.55 0.2 0.4],...
    'Callback',{@doWatershedTransform,wFig}, ...
    'Enable','on','Tag','findEvents');

check_visRect = uicontrol('Style', 'checkbox','Parent',hp_SD,...
    'FontUnits','normalized','Value',1, ...
    'Units','normalized','FontSize',0.75,...
    'String','show detected sparks rect',...
    'Position', [0.115 0.35 0.4 0.1],...
    'Callback',{@showSplittedEvents,wFig});

% slider for sensitivity of watershed transform
txt_spDet = uicontrol('Style','text','FontUnits','normalized',...
    'Parent',hp_SD,'Units','normalized','FontSize',0.7,...
    'Position', [0.025 0.175 0.4 0.1], ...
    'String','watershed sensitivity 50%','Enable','on');

sld_spDet = uicontrol('Style','slider','Parent',hp_SD, ...
    'FontUnits','normalized',...
    'Min',1,'Max',100,'Value',50,'Units','normalized',...
    'Position', [0.025 0.05 0.4 0.1],...
    'Callback', {@sliderFun,wFig},'Enable','on');

uicontrol('Style','text',...
    'Parent',hp_SD,'Units','normalized',...
    'Position', [0.425 0.75 0.225 0.2], ...
    'String',{'min duration of','spark ROI (ms)'},...
    'FontUnits','normalized','FontSize',0.4);
h_edit_MinDurEvent = uicontrol('Style', 'edit','String','20',...
    'Parent',hp_SD,'Units','normalized', ...
    'Position', [0.475 0.575 0.125 0.15],...
    'FontUnits','normalized','FontSize',0.6,...
    'Callback', {@doWatershedTransform,wFig});

uicontrol('Style','text',...
    'Parent',hp_SD,'Units','normalized',...
    'Position', [0.425 0.225 0.225 0.2],...
    'FontUnits','normalized','FontSize',0.4,...
    'String',{'min width of', ['spark (',char(956),'m)']});
h_edit_MinWidthEvent = uicontrol('Style', 'edit','String','1',...
    'Parent',hp_SD,'Units','normalized', ...
    'Position', [0.475 0.05 0.125 0.15],...
    'FontUnits','normalized','FontSize',0.6,...
    'Callback', {@doWatershedTransform,wFig});

% panel to paint mask 
hp_eventMask = uipanel('Title','paint mask','Parent',hp_SD,...
    'Position',[0.625 0 0.15 1],'FontUnits','normalized',...
    'FontSize',0.105,'FontWeight','bold');

% method for spark parameters calculation, create two radio buttons in the button group
eventsMaskPaintTypeRBgroup = uibuttongroup(hp_eventMask, ...
    'Units','normalized',...
    'Title','type of painting:','FontUnits','normalized',...
    'FontSize',0.186,'FontWeight','normal','FontAngle','italic',...
    'Position', [0.025 0.025 0.95 0.475], ...
    'SelectionChangedFcn', '');%  {changeColorOfROI,wFig});

rb1 = uicontrol('Style','radiobutton', ...
    'Parent',eventsMaskPaintTypeRBgroup,...
    'FontUnits','normalized','Value',0,...
    'Units','normalized','FontSize',0.7,...
    'String','add','Position',[0.05 0.55 0.9 0.45]);

rb2 = uicontrol('Style','radiobutton', ...
    'Parent',eventsMaskPaintTypeRBgroup,...
    'FontUnits','normalized','Value',0, ...
    'Units','normalized','FontSize',0.7,...
    'String','remove ','Position',[0.05 0.05 0.9 0.45],'Enable','on');

h_pb_startPaintMask = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> start <br> painting <html>',...
    'FontUnits','normalized','FontSize',0.3, ...
    'Parent',hp_eventMask,'FontWeight','bold',...
    'Units','normalized','Position', [0.025 0.55 0.95 0.4],...
    'Callback',{@startPaintingMask,wFig}, ...
    'Enable','on','Tag','eventsMaskPainting');

% save handles
hObjsW.eventsMaskPaintTypeRBgroup = eventsMaskPaintTypeRBgroup;
hObjsW.rb1AddMask = rb1;
hObjsW.rb2RemoveMask = rb2;


h_push_update = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> update <br> & close <html>',...
    'FontUnits','normalized','FontSize',0.175, ...
    'FontWeight','bold','Parent',hp_SD,...
    'Units','normalized','Position', [0.8 0.2 0.175 0.6],...
    'Callback',{@updateAndClose,wFig},'Enable','on','Tag','findEvents');

hObjsW.check_visRect = check_visRect;
hObjsW.txt_spDet = txt_spDet;
hObjsW.sld_spDet = sld_spDet;

hObjsW.h_edit_MinDurEvent = h_edit_MinDurEvent;
hObjsW.h_edit_MinWidthEvent = h_edit_MinWidthEvent;


%% show event which you want to split and save data

% get data
imgData = getappdata(mainFig, 'imgData');
imgRaw = imgData.imgDataXTfluoR;
imgEvents = imgData.imgDataXTfluoFN;
pxSzX = imgData.pxSzX;
pxSzT = imgData.pxSzT;

% set up axis limits so it is not showing rest of empty mask
% create mask of event
imgEventToSplit = imgData.imgDataXTfluoFN;
bwImg = false(size(imgEventToSplit));
bwImg(eventToSplit.eventDectProp.SubarrayIdx{1}, ...
    eventToSplit.eventDectProp.SubarrayIdx{2}) = true;
imgEventToSplit(~bwImg)=0;

% show in axes
image(imgEventToSplit, ...
    'YData',[0 size(imgEventToSplit,1)], ...
    'XData',[0 size(imgEventToSplit,2)*pxSzT],...
    'CDataMapping','scaled','Parent',ax_img);
set(ax_img,'FontSize',14)
set(get(ax_img,'Ylabel'),'String','x (px)','FontWeight','bold')
set(get(ax_img,'Xlabel'),'String','t (ms)','FontWeight','bold')

% set axis limits
ax_img.XLim = [min(eventToSplit.eventDectProp.SubarrayIdx{2}), ...
    max(eventToSplit.eventDectProp.SubarrayIdx{2})].*pxSzT;
ax_img.YLim = [min(eventToSplit.eventDectProp.SubarrayIdx{1}), ...
    max(eventToSplit.eventDectProp.SubarrayIdx{1})];

% save data
setappdata(wFig, 'mainFig', mainFig)
setappdata(wFig, 'eventToSplit', eventToSplit)
setappdata(wFig, 'hObjsW', hObjsW)
setappdata(wFig, 'imgEventToSplit',imgEventToSplit)


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doWatershedTransform(~,~,wFig)

% get data
hObjsW = getappdata(wFig, 'hObjsW');
mainFig = getappdata(wFig, 'mainFig');
imgData = getappdata(mainFig, 'imgData');
pxSzX = imgData.pxSzX;
pxSzT = imgData.pxSzT;
imgEventToSplit = getappdata(wFig, 'imgEventToSplit');
eventToSplit = getappdata(wFig, 'eventToSplit');

splittedEvents = getappdata(wFig, 'splittedEvents');

% prepare image for watershed
spDetS = round(get(hObjsW.sld_spDet,'Value'));
% parameters of watershed detection/ also filtering, wSp in um and dSp in ms
% filter size is increasing from min values (defined by user) with decreasing sensitivity
wSp = str2double(hObjsW.h_edit_MinWidthEvent.String);
dSp = str2double(hObjsW.h_edit_MinDurEvent.String);

if wSp < pxSzX, wSp = pxSzX; end
if dSp < pxSzT, dSp = pxSzT; end

H = fspecial( 'gaussian', [ceil(wSp/pxSzX), ceil(dSp/pxSzT)], ...
    ceil( min( ceil(wSp/pxSzX), ceil(dSp/pxSzT) )*0.68 ) );
% H = fspecial('average',[ceil(wSp/pxSzX),ceil(dSp/pxSzT)] );
imgEventsWT = imfilter(imgEventToSplit,H,'symmetric');
% create mask of event
bwImg = false(size(imgEventsWT));
bwImg(eventToSplit.eventDectProp.PixelIdxList) = true;
imgEventsWT(~bwImg)=0;
% normalize to [0,1]
imgEventsWT = imgEventsWT ./ max(imgEventsWT(:));

% levels from 1 to 50
imgRange = prctile(imgEventsWT(imgEventsWT~=0), [1 99]);
noiseEst = std(imgEventsWT(imgEventsWT~=0));
nL_min = 1;
nL_max = 5*ceil( (imgRange(2)-imgRange(1)) / noiseEst );
nL = round(((nL_max-nL_min)/(100-1))*spDetS + nL_min - (nL_max-nL_min)/(100-1));
% change image to image with nL levels of intensity
for r = 1:size(imgEventsWT, 1)
    for c = 1:size(imgEventsWT, 2)
        if imgEventsWT(r,c) > 0
            imgEventsWT(r,c) = round(imgEventsWT(r,c)*nL);
        end
    end
end

% watershed
L = watershed(-imgEventsWT,8);
L(L>0)=1;
bwImg(~logical(L)) = 0;

% figure
% subplot(2,1,1)
% imagesc(imgEventsWT)
% subplot(2,1,2)
% imagesc(bwImg)


% merge with manualy selected mask if exist
if isfield(getappdata(wFig),'maskPaint')
    maskPaint = getappdata(wFig,'maskPaint');

    % find difference
    diffMasks = xor(bwImg,maskPaint.tempMaskOfEvents);
    % find additions
    maskAdd = maskPaint.tempMaskOfEvents==true;
    maskAdd = maskAdd & diffMasks;
    % find removal
    maskRemoval = maskPaint.tempMaskOfEvents==false;
    maskRemoval = maskRemoval & diffMasks;
    
    % create final mask
    bwImg = bwImg | maskAdd;
    bwImg = bwImg & ~maskRemoval;

end

% find all connected components
CC_event = bwconncomp(bwImg,8);

% calculate properties of all CC (regions)
statEvents = regionprops(CC_event,imgEventToSplit, 'Centroid','PixelList',...
    'PixelIdxList','SubarrayIdx','PixelValues','MaxIntensity',...
    'WeightedCentroid','BoundingBox','Image','ConvexHull','Area');
% calculate boundaries of events
eventsBoundaries = bwboundaries(bwImg,8,'noholes');
% add them to stat events
[statEvents.('Boundary')] = eventsBoundaries{:};

% show newly detected events after spliting original event

% remove previous
if ~isempty(splittedEvents)
    if isfield(splittedEvents,'detectedEventsRec')
        delete(findall(hObjsW.ax_img,'Tag','detectedEventRecText'))
        delete(splittedEvents.detectedEventsRec)
        delete(splittedEvents.detectedEventsMask)
        rmappdata(wFig,'splittedEvents')
    end
end

% create uicontextmenu to change type of event later
% 1 = spark, 2 = long-lasting spark, 3 = mini-wave, 4 = wave, 5 = transient
cm = uicontextmenu(wFig);
m1 = uimenu(cm,'Text','1. spark', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEventLocal, wFig});
m2 = uimenu(cm,'Text','2. LL spark', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEventLocal, wFig});
m3 = uimenu(cm,'Text','3. mini-wave', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEventLocal, wFig});
m4 = uimenu(cm,'Text','4. wave', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEventLocal, wFig});
m5 = uimenu(cm,'Text','5. transient', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEventLocal, wFig});
m6 = uimenu(cm,'Text','6. caffeine transient', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEventLocal, wFig});

detectedEventsRec = gobjects(length(statEvents),1);
detectedEventsMask = gobjects(length(statEvents),1);
maskOfAcceptedSparks = true(1,length(statEvents));
% plot bounding rectangles
for i = 1:length(statEvents)
    
    pos = statEvents(i).BoundingBox;
    % check dimensions of events
    if ( pos(3)*pxSzT < dSp ) || ( pos(4)*pxSzX < dSp )
        eventColor = [0 0 0];
        maskOfAcceptedSparks(i) = false;
    else
        eventColor = eventToSplit.c;
    end
    
    detectedEventsRec(i) = rectangle( ...
        'Position',[pos(1)*pxSzT pos(2) pos(3)*pxSzT pos(4)],...
        'Parent',hObjsW.ax_img,'EdgeColor',eventColor,...
        'LineWidth',2,'Tag',num2str(i),...
        'ButtonDownFcn', {@setTypeOfCalciumEventLocal, wFig}, ...
        'UIContextMenu', cm);
    
    % show detected mask
    detectedEventsMask(i) = patch( ...
        'Faces',(1:numel(statEvents(i).Boundary(:,2))), ...
        'Vertices',[statEvents(i).Boundary(:,2).*pxSzT, ...
                    statEvents(i).Boundary(:,1)], ...
        'FaceColor',eventColor, 'FaceAlpha',0.2, ...
        'EdgeColor',eventColor, 'EdgeAlpha',1, ...
        'Parent',hObjsW.ax_img, ...
        'LineWidth',1, 'Tag',num2str(i),...
        'ButtonDownFcn', {@setTypeOfCalciumEventLocal, wFig}, ...
        'UIContextMenu', cm);
    
    text(hObjsW.ax_img,detectedEventsRec(i).Position(1),...
        detectedEventsRec(i).Position(2),num2str(i),...
        'FontSize',20,'VerticalAlignment','bottom','FontWeight','bold',...
        'Color',eventColor,'Tag','detectedEventRecText')
    
    clearvars pos
end

splittedEvents.detectedEvents = statEvents;
splittedEvents.detectedEventsRec = detectedEventsRec;
splittedEvents.detectedEventsMask = detectedEventsMask;
splittedEvents.maskOfAcceptedSparks = maskOfAcceptedSparks;
% save also type of detected events, used later when creating dataset
% for machine learning
splittedEvents.typeOfEvent = ones(1,length(statEvents)) .* eventToSplit.type;

% save data
setappdata(wFig,'splittedEvents',splittedEvents);
set(hObjsW.check_visRect,'Value',1)

end


function updateAndClose(~,~,wFig)
% update mask and detected events structure in main figure
% get data
eventToSplit = getappdata(wFig, 'eventToSplit');
splittedEvents = getappdata(wFig, 'splittedEvents');

mainFig = getappdata(wFig, 'mainFig');
sparkDetection = getappdata(mainFig, 'sparkDetection');
hObjs = getappdata(mainFig, 'hObjs');
imgData = getappdata(mainFig, 'imgData');
pxSzT = imgData.pxSzT;

% % delete original events, and redraw them with splitted one
% delete(sparkDetection.detectedEventsRec)
% delete(sparkDetection.detectedEventsMask)
% delete(findall(hObjs.ax_img, 'Tag','detectedEventRecText'))

% update sparkDetection structure with splitted events
% do not take rejected events
acceptedEvents = structfun( ...
    @(x) x(splittedEvents.maskOfAcceptedSparks), splittedEvents, ...
    'UniformOutput', 0 );

newDetectedEvents = [ ...
    sparkDetection.detectedEvents(1:eventToSplit.idx-1); ...
    acceptedEvents.detectedEvents; ...
    sparkDetection.detectedEvents(eventToSplit.idx+1:end) ];
newTypeOfEvent = [ ...
    sparkDetection.typeOfEvent(1:eventToSplit.idx-1), ...
    acceptedEvents.typeOfEvent, ...
    sparkDetection.typeOfEvent(eventToSplit.idx+1:end) ];
newMaskOfAcceptedSparks = [ ...
    sparkDetection.maskOfAcceptedSparks(1:eventToSplit.idx-1), ...
    acceptedEvents.maskOfAcceptedSparks, ...
    sparkDetection.maskOfAcceptedSparks(eventToSplit.idx+1:end) ];

sparkDetection.detectedEvents = newDetectedEvents;
sparkDetection.typeOfEvent = newTypeOfEvent;
sparkDetection.maskOfAcceptedSparks = newMaskOfAcceptedSparks;

if isfield(sparkDetection, ...
        'analyzedEvntsBrowserTbl')
    % remove it, and reanalyze events after
    % you are done with deleting/splitting events
    close(findall(0, 'Type','Figure', ...
        'Tag','CaEventsBrowser'))
    sparkDetection = rmfield(sparkDetection, ...
        "analyzedEvntsBrowserTbl");
end

% save changes
setappdata(mainFig, 'sparkDetection', sparkDetection)
% update detected events
eventsDetection([], [], mainFig, 'update')

% 
% % create uicontextmenu to change type of event later
% % 1 = spark, 2 = long-lasting spark, 3 = mini-wave, 4 = wave, 5 = transient
% cm = uicontextmenu(mainFig);
% m1 = uimenu(cm,'Text','1. spark', ...
%     'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
% m2 = uimenu(cm,'Text','2. LL spark', ...
%     'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
% m3 = uimenu(cm,'Text','3. mini-wave', ...
%     'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
% m4 = uimenu(cm,'Text','4. wave', ...
%     'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
% m5 = uimenu(cm,'Text','5. transient', ...
%     'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
% m6 = uimenu(cm,'Text','6. caffeine transient', ...
%     'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
% m7 = uimenu(cm,'Text','split event (watershed)', ...
%         'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
% m8 = uimenu(cm,'Text','show parameters of detected event', ...
%         'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
% 
% % show newly splitted in mainFig axes
% detectedEventsRec = gobjects(length(newDetectedEvents),1);
% detectedEventsMask = gobjects(length(newDetectedEvents),1);
% % plot bounding rectangles
% for i = 1:length(newDetectedEvents)
% 
%     pos = newDetectedEvents(i).BoundingBox;
%     % check dimensions of events
%     if newMaskOfAcceptedSparks(i)
%         switch newTypeOfEvent(i)
%             case 1
%                 eventColor = [1 0 0];
%             case 2
%                 eventColor = [1 0 1];
%             case 3
%                 eventColor = [0 1 1];
%             case 4
%                 eventColor = [0 1 0];
%             case 5
%                 eventColor = [1 1 1];
%             case 6
%                 eventColor = [1 1 0];
%         end
%     else
%         eventColor = [0 0 0];
%     end
% 
%     detectedEventsRec(i) = rectangle( ...
%         'Position',[pos(1)*pxSzT pos(2) pos(3)*pxSzT pos(4)],...
%         'Parent',hObjs.ax_img,'EdgeColor',eventColor,...
%         'LineWidth',3,'Tag',num2str(i),...
%         'ButtonDownFcn', {@setTypeOfCalciumEvent, mainFig}, ...
%         'UIContextMenu', cm);
% 
%     % show detected mask
%     detectedEventsMask(i) = patch( ...
%         'Faces',(1:numel(newDetectedEvents(i).Boundary(:,2))), ...
%         'Vertices',[newDetectedEvents(i).Boundary(:,2).*pxSzT, ...
%                     newDetectedEvents(i).Boundary(:,1)], ...
%         'FaceColor',eventColor, 'FaceAlpha',0.05, ...
%         'EdgeColor',eventColor, 'EdgeAlpha',1, ...
%         'Parent',hObjs.ax_img, ...
%         'LineWidth',2, 'Tag',num2str(i),...
%         'ButtonDownFcn', {@setTypeOfCalciumEvent, mainFig}, ...
%         'UIContextMenu', cm);
% 
%     text(hObjs.ax_img,detectedEventsRec(i).Position(1),...
%         detectedEventsRec(i).Position(2),num2str(i),...
%         'FontSize',10,'VerticalAlignment','bottom','FontWeight','bold',...
%         'Color',eventColor,'Tag','detectedEventRecText')
% 
%     clearvars pos
% end
% 
% % update data in mainFig
% sparkDetection.detectedEvents = newDetectedEvents;
% sparkDetection.typeOfEvent = newTypeOfEvent;
% sparkDetection.maskOfAcceptedSparks = newMaskOfAcceptedSparks;
% sparkDetection.detectedEventsRec = detectedEventsRec;
% sparkDetection.detectedEventsMask = detectedEventsMask;
% % recalculate spark frequency
% calcSparkFreq(mainFig, 1);
% sparkDetection.sparkFreq = str2double(hObjs.txt_correctedSpF.String{1});
% setappdata(mainFig,'sparkDetection',sparkDetection)

% close window
delete(wFig)

end


function sliderFun(hObj,~,wFig)

% update text field of slider
hObjsW = getappdata(wFig, 'hObjsW');

hObjsW.txt_spDet.String = ...
    sprintf('watershed sensitivity %d %%', round(hObj.Value));

% do watershed
doWatershedTransform([],[],wFig)

end


function showSplittedEvents(hObj,~,wFig)

val = get(hObj,'Value');
splittedEvents = getappdata(wFig,'splittedEvents');
hObjsW = getappdata(wFig,'hObjsW');
rectTxt = findall(hObjsW.ax_img,'Tag','detectedEventRecText');

if isfield(splittedEvents,'detectedEventsRec')
    
    if val == 1
        set(splittedEvents.detectedEventsRec,'Visible','on')
        set(splittedEvents.detectedEventsMask,'Visible','on')
        set(rectTxt,'Visible','on')
    else
        set(splittedEvents.detectedEventsRec,'Visible','off')
        set(splittedEvents.detectedEventsMask,'Visible','off')
        set(rectTxt,'Visible','off')
    end
    
else
    set(hObj,'Value',0)
    return   
end

end


function setTypeOfCalciumEventLocal(hObj,E,wFig)
%%%%%%%%%%%%%
% with click of mouse button 1 change state of detected event between
% accepted (color) and rejected (black rectangle)
%%%%%%%%%%%%%
% with click of mouse button 2, call contextmenu and change state of detected event among type of
% event: 1 = spark (red), 
% 2 = long-lasting spark (magenta), 
% 3 = mini-wave (cyan), 
% 4 = wave (green), 
% 5 = transient (white)
% 6 = caffeine transient (yellow)

% get data
hObjsW = getappdata(wFig, 'hObjsW');
splittedEvents = getappdata(wFig,'splittedEvents');

% handle interactions
switch hObj.Type
    
    case {'rectangle', 'patch'}
        idx = str2double(hObj.Tag);
        typeOfEvent = splittedEvents.typeOfEvent(idx);
        h_rectTxt = findobj(hObjsW.ax_img,'Type','text', ...
                'Tag','detectedEventRecText', 'String', hObj.Tag);
        
        switch hObj.Type
            case 'rectangle'
                h_rect = hObj; 
                h_patch = findobj(hObjsW.ax_img,'Type','patch', ...
                    'Tag',hObj.Tag);
                
            case 'patch'
                h_patch = hObj;
                h_rect = findobj(hObjsW.ax_img,'Type','rectangle', ...
                    'Tag',hObj.Tag);
     
        end
               
        switch E.Button
            
            case 1
                if any(hObj.EdgeColor ~= [0 0 0])
                    state = false;
                    c = [0 0 0];
                    
                else
                    state = true;
                    switch typeOfEvent
                        case 1
                            c = [1 0 0];
                        case 2
                            c = [1 0 1];
                        case 3
                            c = [0 1 1];
                        case 4
                            c = [0 1 0];
                        case 5
                            c = [1 1 1];
                        case 6
                            c = [1 1 0];
                    end

                end
        end
        
    case 'uimenu'
        
        % handle to current object
        hCurrentObj = gco(wFig);
        switch hCurrentObj.Type
            case 'rectangle'
                h_rect = hCurrentObj; 
                h_patch = findobj(hObjsW.ax_img,'Type','patch', ...
                    'Tag',hCurrentObj.Tag);
                
            case 'patch'
                h_patch = hCurrentObj;
                h_rect = findobj(hObjsW.ax_img,'Type','rectangle', ...
                    'Tag',hCurrentObj.Tag);
     
        end
  
        idx = str2double(h_rect.Tag);
        state = true;
        h_rectTxt = findobj(hObjsW.ax_img,'Type','text', ...
                'Tag','detectedEventRecText', 'String', h_rect.Tag);

        switch hObj.Text
            
            case '1. spark'
                typeOfEvent = 1;
                c = [1 0 0];
            
            case '2. LL spark'
                typeOfEvent = 2;
                c = [1 0 1];
                
            case '3. mini-wave'
                typeOfEvent = 3;
                c = [0 1 1];
                
            case '4. wave'
                typeOfEvent = 4;
                c = [0 1 0];
                
            case '5. transient'
                typeOfEvent = 5;
                c = [1 1 1];
                
            case '6. caffeine transient'
                typeOfEvent = 6;
                c = [1 1 0];
                       
        end
    
end

% change properties
try
    % set up rectangle
    h_rect.EdgeColor = c;
    h_rectTxt.Color = c;
    % set up mask
    h_patch.EdgeColor = c;
    h_patch.FaceColor = c;

    splittedEvents.maskOfAcceptedSparks(idx) = state;
    splittedEvents.typeOfEvent(idx) = typeOfEvent;

    % save changes
    setappdata(wFig, 'splittedEvents', splittedEvents)
catch 
   return
end    
    
end


function startPaintingMask(hObj,E,wFig)
% when pressed create rectangle and start painting
hObjsW = getappdata(wFig,'hObjsW');
splittedEvents = getappdata(wFig,'splittedEvents');
imgEventToSplit = getappdata(wFig,'imgEventToSplit');
eventToSplit = getappdata(wFig, 'eventToSplit');
mainFig = getappdata(wFig, 'mainFig');
imgData = getappdata(mainFig, 'imgData');
pxSzT = imgData.pxSzT;

switch hObj.String
    case '<html> <p align="center"> start <br> painting <html>'
        
        % change name
        hObj.String = '<html> <p align="center"> stop <br> painting <html>';
 
        % delete rectangle
        if ~isempty(splittedEvents)
            delete(splittedEvents.detectedEventsMask)
            delete(splittedEvents.detectedEventsRec)
            delete(findall(hObjsW.ax_img,'Tag','detectedEventRecText'))
            % create mask of events
            tempMaskOfEvents = false(size(imgEventToSplit));
            for i = 1:length(splittedEvents.detectedEvents)
                
                tempMaskOfEvents( ...
                    splittedEvents.detectedEvents(i).SubarrayIdx{1}, ...
                    splittedEvents.detectedEvents(i).SubarrayIdx{2} ) = ...
                    splittedEvents.detectedEvents(i).Image ;
            end
            
        else
            % create mask of event
            tempMaskOfEvents = false(size(imgEventToSplit));
            tempMaskOfEvents(eventToSplit.eventDectProp.PixelIdxList) = true;
            
        end
        
        % create dummy axes only for selecting mask
        ax_maskROI = axes('Parent',wFig,...
            'Position',hObjsW.ax_img.Position);
        set(ax_maskROI,'XTick',[],'YTick',[],'Color','none','Tag','tempROIaxes')
        % disable zooming, moving by mouse in axes
        disableDefaultInteractivity(ax_maskROI)
        % disableDefaultInteractivity(hObjsA.ax_detEvent) % maybe
        
        % create temporal transparent mask
        h_tempMaskImg = image(ax_maskROI, ...
            'CData',tempMaskOfEvents,...
            'YData',[1 size(tempMaskOfEvents,1)], ...
            'XData',[0 size(tempMaskOfEvents,2)*pxSzT],...
            'CDataMapping','scaled','AlphaData',0.25);
        ax_maskROI.YAxis.Direction = 'reverse';
        
        % set axis limits
        ax_maskROI.XLim = [min(eventToSplit.eventDectProp.SubarrayIdx{2}), ...
            max(eventToSplit.eventDectProp.SubarrayIdx{2})].*pxSzT;
        ax_maskROI.YLim = [min(eventToSplit.eventDectProp.SubarrayIdx{1}), ...
            max(eventToSplit.eventDectProp.SubarrayIdx{1})];
        
        % create ROI to paint mask
        hROI = drawrectangle(ax_maskROI);
        switch hObjsW.eventsMaskPaintTypeRBgroup.SelectedObject.String
            case 'add'
                hROI.Color = [0 1 0];
            case 'remove '
                hROI.Color = [1 0 0];
        end
        
        % notify during interaction with ROI
        notify(hROI,'MovingROI')
        
        % listeners/callbacks % 'MovingROI' ROIMoved
        lh_moving = addlistener(hROI,'MovingROI', ...
            @updateMaskOfEventWithROI);
        
        % save data
        maskPaint.ax_maskROI = ax_maskROI;
        maskPaint.tempMaskOfEvents = tempMaskOfEvents;
        maskPaint.h_tempMaskImg = h_tempMaskImg;
        maskPaint.hROI = hROI;
        setappdata(wFig,'maskPaint',maskPaint)

    case '<html> <p align="center"> stop <br> painting <html>'
        
        maskPaint = getappdata(wFig,'maskPaint');
        %eventToSplit = getappdata(wFig, 'eventToSplit');
        % delete ROI
        delete(maskPaint.hROI)
        % delete axes for temporal mask
        delete(maskPaint.ax_maskROI)
        % do watershed
        doWatershedTransform([],[],wFig)

        % change name
        hObj.String = '<html> <p align="center"> start <br> painting <html>';

end

end



function updateMaskOfEventWithROI(hRect,E)

wFig = hRect.Parent.Parent;
hObjsW = getappdata(wFig,'hObjsW');
maskPaint = getappdata(wFig,'maskPaint');
h_tempMaskImg = maskPaint.h_tempMaskImg;
tempMask = maskPaint.tempMaskOfEvents;
  
% get mask of rectangle
rectMask = hRect.createMask(h_tempMaskImg);

% update temp mask
switch hObjsW.eventsMaskPaintTypeRBgroup.SelectedObject.String
    case 'add'
        tempMask = tempMask | rectMask;
    case 'remove '
        tempMask = tempMask & (~rectMask);
end
 
% update image
h_tempMaskImg.CData = tempMask;

% save
maskPaint.tempMaskOfEvents = tempMask;
setappdata(wFig,'maskPaint',maskPaint)

end



function varargout = changeColorOfROI(varargin)
keyboard
%hObjsW = getappdata(wFig,'hObjsW');

keyboard

% create ROI to paint mask
hROI = drawrectangle(ax_maskROI);
switch hObjsW.eventsMaskPaintTypeRBgroup.SelectedObject.String
    case 'add'
        hROI.Color = [0 1 0];
    case 'remove '
        hROI.Color = [1 0 0];
end

end


