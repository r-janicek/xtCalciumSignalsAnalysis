function eventsDetection(h_O, ~, mainFig, detectionFlag)
% detect calcium events in linescan using algorithm developed by:
% Detecting Ca2+ sparks on stationary and varying baselines
% Peter Bankhead, C. Norman Scholfield, Tim M. Curtis, and J. Graham McGeown
% American Journal of Physiology-Cell Physiology 2011 301:3, C717-C728

% detectionFlag = 'new' or 'update' detection 
% (update is calling when event is deleted or splitted)

set(mainFig,'Pointer','watch')
drawnow
% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
imgDataXTfluoR = imgData.imgDataXTfluoR;
%imgDataXTfluoRN = imgData.imgDataXTfluoRN;
pxSzT = imgData.pxSzT;
sparkDetection = getappdata(mainFig,'sparkDetection');

% remove previous
if ~isempty(sparkDetection)
    if isfield(sparkDetection,'detectedEventsRec')
        delete(findall(hObjs.ax_img,'Tag','detectedEventRecText'))
        delete(sparkDetection.detectedEventsRec)
        delete(sparkDetection.detectedEventsMask)
        rmappdata(mainFig,'sparkDetection')
    end
end
% get events
switch detectionFlag
    case 'new'
        % last pressed pushbutton to detect events
        setappdata(mainFig,'lastPressedPusbutton',h_O);
        % find events
        events = findEvents(mainFig, imgDataXTfluoR);
    case 'update'
        events = sparkDetection.detectedEvents;
end
% create uicontextmenu to change type of event later
% 1 = spark, 2 = long-lasting spark, 3 = mini-wave, 4 = wave, 5 = transient
cm = uicontextmenu(mainFig);
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
m7 = uimenu(cm,'Text','split event (watershed)', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
m8 = uimenu(cm,'Text','show parameters of detected event', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
m9 = uimenu(cm,'Text','delete detected event', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
m10 = uimenu(cm,'Text','delete multiple detected events', ...
    'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});


detectedEventsRec = gobjects(length(events),1);
detectedEventsMask = gobjects(length(events),1);
% plot bounding rectangles
for i = 1:length(events)
    pos = events(i).BoundingBox;
    if ~isempty(sparkDetection)
        if isfield(sparkDetection,'maskOfAcceptedSparks')
            if sparkDetection.maskOfAcceptedSparks(i)
                switch sparkDetection.typeOfEvent(i)
                    case 1
                        eventColor = [1 0 0];
                    case 2
                        eventColor = [1 0 1];
                    case 3
                        eventColor = [0 1 1];
                    case 4
                        eventColor = [0 1 0];
                    case 5
                        eventColor = [1 1 1];
                    case 6
                        eventColor = [1 1 0];
                end
            else
                eventColor = [0 0 0];
            end
        end
    else
        eventColor = [1 0 0];
    end
    detectedEventsRec(i) = rectangle( ...
        'Position',[pos(1)*pxSzT pos(2) pos(3)*pxSzT pos(4)],...
        'Parent',hObjs.ax_img,'EdgeColor',eventColor,...
        'LineWidth',2,'Tag',num2str(i),...
        'ButtonDownFcn', {@setTypeOfCalciumEvent, mainFig}, ...
        'UIContextMenu', cm);

    % show detected mask
    detectedEventsMask(i) = patch( ...
        'Faces',(1:numel(events(i).Boundary(:,2))), ...
        'Vertices',[events(i).Boundary(:,2).*pxSzT, ...
        events(i).Boundary(:,1)], ...
        'FaceColor',eventColor, 'FaceAlpha',0.05, ...
        'EdgeColor',eventColor, 'EdgeAlpha',1, ...
        'Parent',hObjs.ax_img, ...
        'LineWidth',1, 'Tag',num2str(i),...
        'ButtonDownFcn', {@setTypeOfCalciumEvent, mainFig}, ...
        'UIContextMenu', cm);

    text(hObjs.ax_img,detectedEventsRec(i).Position(1),...
        detectedEventsRec(i).Position(2),num2str(i),...
        'FontSize',10,'VerticalAlignment','bottom','FontWeight','bold',...
        'Color',detectedEventsRec(i).EdgeColor,'Tag','detectedEventRecText')

    clearvars pos
end

sparkDetection.detectedEvents = events;
sparkDetection.detectedEventsRec = detectedEventsRec;
sparkDetection.detectedEventsMask = detectedEventsMask;
% save also type of detected events, used later when creating dataset
% for machine learning
switch detectionFlag
    case 'new'
        sparkDetection.typeOfEvent = ones(1,length(events));
        sparkDetection.maskOfAcceptedSparks  = true(1,length(events));
    case 'update'

end

setappdata(mainFig,'sparkDetection',sparkDetection);
set(hObjs.check_visRect, 'Value',1)

% calculate frequency of events in image
calcSparkFreq(mainFig, false)

set(mainFig,'Pointer','arrow')
drawnow

% if analysis type is to detect sparks then do also image normalization
switch getappdata(mainFig, 'analysisType')
    case 'spark detection'
        switch detectionFlag
            case 'new'
                norm_calc([], [], mainFig, [])
            case 'update'
        end
    otherwise    
end
end

