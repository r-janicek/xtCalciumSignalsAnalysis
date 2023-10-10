function setCriteriaSparkDetection(hO,~,mainFig)

hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');
pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;

% get data; get parameters and check theirs values
FWHM = str2double(get(hObjs.h_edit_fFWHM,'String'));
smoothIter = str2double(get(hObjs.h_edit_smoothIter,'String'));
baseIter = str2double(get(hObjs.h_edit_baseIter,'String'));
tresh = str2double(get(hObjs.h_edit_tresh,'String'));
expFactor = str2double(get(hObjs.h_edit_expFactor,'String'));

if baseIter > 10 
    baseIter = 10;
    set(hObjs.h_edit_baseIter,'String',num2str(baseIter))
end

if smoothIter >= baseIter
    smoothIter = baseIter-1;
    set(hObjs.h_edit_smoothIter,'String',num2str(smoothIter))
end
      
profileAnalysis = getappdata(mainFig,'profileAnalysis');
% last pressed pushbutton to detect events
detectEventsPB = getappdata(mainFig,'lastPressedPusbutton'); 

if strcmp(detectEventsPB.Tag,'findEvents')
    eventsDetection(hObjs.h_push_sparks,[],mainFig,'new')
else
    
    cropRoi = profileAnalysis.croppedDataWholeAreaRoi; % filtred and possibly normalized data
    cropRoiR = profileAnalysis.croppedDataWholeAreaRoiR; % raw data
    
    imgRaw = profileAnalysis.croppedDataWholeAreaRoiR;   
    statEventsNew = findEvents(mainFig,cropRoiR);
    
    % plot new crosses and boundaries rectangles of detected events
    % delete previous 
    statEventsOld = profileAnalysis.statEventsSpRec;
    delete([statEventsOld.centreLine])
    delete([statEventsOld.eventRec])
    
    % find centre of the event
    pos_centre_x = arrayfun(@(x) x.WeightedCentroid(1,1),statEventsNew,'UniformOutput',1);
    pos_centre_y = round(arrayfun(@(x) x.WeightedCentroid(1,2),statEventsNew,'UniformOutput',1));
    
    %pos_max = num2cell([pos_centre_y,pos_centre_x],2);
    
    for i = 1:length(statEventsNew)
        
        event = cropRoi(statEventsNew(i).PixelIdxList);
        [~,p_m] = max(event);
        px = statEventsNew(i).PixelIdxList;
        
        col(i,1) = ceil(px(p_m)/size(cropRoi,1));
        row = rem(px(p_m),size(cropRoi,1));
        
    end
    
    try
        pos_max = num2cell([pos_centre_y,col],2);
    catch
        pos_max = [];
    end
    
    % show centre and boundary rectangle of detected events
    for i=1:length(statEventsNew)
        
        % bounding rectangle
        pos = statEventsNew(i).BoundingBox;
        eventRec = rectangle('Position',[pos(1)*pxSzT pos(2) pos(3)*pxSzT pos(4)],...
            'Parent',hObjs.ax_img_sparks,'EdgeColor','r','LineStyle',':','LineWidth',2);
        
        % centre line
        centre = pos_max{i};
        centreLine = line([centre(2)*pxSzT-pxSzT centre(2)*pxSzT-pxSzT],[centre(1) centre(1)],...
            'Parent',hObjs.ax_img_sparks,'Color','r','LineStyle','none',...
            'Marker','+','MarkerSize',20,'LineWidth',3,'PickableParts','all',...
            'ButtonDownFcn',{@profileROIButtonDownFcn,mainFig});
        
        [statEventsNew(i).centreLine] = centreLine;
        [statEventsNew(i).eventRec] = eventRec;
        
    end
        
    profileAnalysis.statEventsSpRec = statEventsNew;
    setappdata(mainFig,'profileAnalysis',profileAnalysis)

end
    

end

