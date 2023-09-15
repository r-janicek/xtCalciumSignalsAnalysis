function eventsDetection(h_O,~,mainFig)

set(mainFig,'Pointer','watch')
drawnow

% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
imgDataXTfluoR = imgData.imgDataXTfluoR;
%imgDataXTfluoRN = imgData.imgDataXTfluoRN;

pxSzT = imgData.pxSzT;

sparkDetection = getappdata(mainFig,'sparkDetection');
visRect = get(hObjs.check_visRect,'Value');

% if strcmp(get(h_O,'Tag'),'findEvents')
       
    % remove previous
    if ~isempty(sparkDetection)
        if isfield(sparkDetection,'detectedEventsRec') 
            delete(findall(hObjs.ax_img,'Tag','detectedEventRecText'))
            delete(sparkDetection.detectedEventsRec)
            delete(sparkDetection.detectedEventsMask)
            rmappdata(mainFig,'sparkDetection')           
        end       
    end
    
    % last pressed pushbutton to detect events
    setappdata(mainFig,'lastPressedPusbutton',h_O);
    
    % find events  
    events = findEvents(mainFig, imgDataXTfluoR);

    % create uicontextmenu to change type of event later
    % 1 = spark, 2 = long-lasting spark, 3 = mini-wave, 4 = wave, 5 = transient
    cm = uicontextmenu(mainFig);
    m1 = uimenu(cm,'Text','1. spark', ...
        'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
    m2 = uimenu(cm,'Text','2. LL spark', ...
        'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
    m3 = uimenu(cm,'Text','3. mini-wave', ...
        'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
    m4 = uimenu(cm,'Text','4. wave', ...
        'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
    m5 = uimenu(cm,'Text','5. transient', ...
        'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
    m6 = uimenu(cm,'Text','6. caffeine transient', ...
        'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
    m7 = uimenu(cm,'Text','split event (watershed)', ...
        'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
    m8 = uimenu(cm,'Text','show parameters of detected event', ...
        'MenuSelectedFcn',{@setTypeOfCalciumEvent, mainFig});
    
    detectedEventsRec = gobjects(length(events),1);
    detectedEventsMask = gobjects(length(events),1);
    % plot bounding rectangles
    for i = 1:length(events)
 
        pos = events(i).BoundingBox;
        detectedEventsRec(i) = rectangle( ...
            'Position',[pos(1)*pxSzT pos(2) pos(3)*pxSzT pos(4)],...
            'Parent',hObjs.ax_img,'EdgeColor','r',...
            'LineWidth',2,'Tag',num2str(i),...
            'ButtonDownFcn', {@setTypeOfCalciumEvent, mainFig}, ...
            'UIContextMenu', cm);
       
        % show detected mask
        detectedEventsMask(i) = patch( ...
            'Faces',(1:numel(events(i).Boundary(:,2))), ...
            'Vertices',[events(i).Boundary(:,2).*pxSzT, ...
                        events(i).Boundary(:,1)], ...
            'FaceColor',[1 0 0], 'FaceAlpha',0.05, ...
            'EdgeColor',[1 0 0], 'EdgeAlpha',1, ...
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
    sparkDetection.typeOfEvent = ones(1,length(events));
    sparkDetection.maskOfAcceptedSparks  = true(1,length(events));
    
    setappdata(mainFig,'sparkDetection',sparkDetection);
    set(hObjs.check_visRect,'Value',1)
       
    % calculate frequency of events in image
    calcSparkFreq(mainFig,false)    
                 
% else    
%     
%     keyboard
%     
%     if ~isfield(getappdata(mainFig),'S_event')
%         set(getappdata(mainFig,'main_fig'),'Pointer','arrow')
%         drawnow
%         return
%     end
%  
%     events = getappdata(mainFig,'S_event');  
%     eventParams = calcEventParams(Img_data,events,[],[],[],mainFig);               
%     setappdata(mainFig,'calciumSparksParameters',eventParams);
%        
% end

set(mainFig,'Pointer','arrow')
drawnow

% if analysis type is to detect sparks then do also image normalization
switch getappdata(mainFig,'analysisType')
    case 'spark detection'
        norm_calc_sparksDet(mainFig)
        
    otherwise    
end


end

