function setTypeOfCalciumEvent(hObj,E,mainFig)
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
hObjs = getappdata(mainFig, 'hObjs');
detectedEvnts = getappdata(mainFig,'sparkDetection');

% handle interactions
switch hObj.Type
    
    case {'rectangle', 'patch'}
        idx = str2double(hObj.Tag);
        typeOfEvent = detectedEvnts.typeOfEvent(idx);
        h_rectTxt = findobj(hObjs.ax_img,'Type','text', ...
                'Tag','detectedEventRecText', 'String', hObj.Tag);
        
        switch hObj.Type
            case 'rectangle'
                h_rect = hObj; 
                h_patch = findobj(hObjs.ax_img,'Type','patch', ...
                    'Tag',hObj.Tag);
                
            case 'patch'
                h_patch = hObj;
                h_rect = findobj(hObjs.ax_img,'Type','rectangle', ...
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
        hCurrentObj = gco(mainFig);
        switch hCurrentObj.Type
            case 'rectangle'
                h_rect = hCurrentObj; 
                h_patch = findobj(hObjs.ax_img,'Type','patch', ...
                    'Tag',hCurrentObj.Tag);
                
            case 'patch'
                h_patch = hCurrentObj;
                h_rect = findobj(hObjs.ax_img,'Type','rectangle', ...
                    'Tag',hCurrentObj.Tag);
        end
  
        idx = str2double(h_rect.Tag);
        state = true;
        h_rectTxt = findobj(hObjs.ax_img,'Type','text', ...
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
            
            case 'delete detected event'
                % answer = questdlg( ...
                %     'Would you like to delete selected event?', ...
                % 	'Delete selected event.', ...
                % 	'YES','NO','NO');
                answer = 'YES';
                switch answer
                    case 'YES'
                        % delete event from detection output
                        newDetectedEvents = detectedEvnts.detectedEvents;
                        newDetectedEvents(idx) = [];
                        newTypeOfEvent = detectedEvnts.typeOfEvent;
                        newTypeOfEvent(idx) = [];
                        newMaskOfAcceptedSparks = detectedEvnts.maskOfAcceptedSparks;
                        newMaskOfAcceptedSparks(idx) = [];
                        if isfield(detectedEvnts, ...
                                'analyzedEvntsBrowserTbl')
                            % remove it, and reanalyze events after 
                            % you are done with deleting/splitting events 
                            close(findall(0, 'Type','Figure', ...
                                'Tag','CaEventsBrowser'))
                            detectedEvnts = rmfield(detectedEvnts, ...
                                "analyzedEvntsBrowserTbl");
                        end
                        detectedEvnts.detectedEvents = newDetectedEvents;
                        detectedEvnts.typeOfEvent = newTypeOfEvent;
                        detectedEvnts.maskOfAcceptedSparks = newMaskOfAcceptedSparks;
                        % save changes
                        setappdata(mainFig, 'sparkDetection', detectedEvnts)
                        % update detected events
                        eventsDetection(hObj, [], mainFig, 'update')

                    case 'NO'
                        return
                end


            case 'delete multiple detected events'
                h_img = findobj(hObjs.ax_img, 'Type', 'Image');
                h_img.PickableParts = 'none';
                % change pointer
                mainFig.Pointer = 'cross';
                drawnow
                hObjs.ax_img.ButtonDownFcn = ...
                    {@deleteMultipleDetectedEvents};
                mainFig.Pointer = 'arrow';
                drawnow
                

            case 'split event (watershed)'
                % apply watershed transform on selected event only
                % get data of event
                eventToSplit = struct( ...
                    'idx', idx, ...
                    'eventDectProp', detectedEvnts.detectedEvents(idx), ...
                    'type', detectedEvnts.typeOfEvent(idx), ...
                    'state', detectedEvnts.maskOfAcceptedSparks(idx), ...
                    'c', h_rect.EdgeColor, ...
                    'h_rect' , h_rect, ...
                    'h_patch', h_patch);
                splitSelectedEventWindow(mainFig, eventToSplit)


            case 'show parameters of detected event'
                imgData = getappdata(mainFig,'imgData');
                % selected method of parameters calculation
                switch hObjs.calcSpParamMethodRBgroup.SelectedObject.String
                    case '2D gaussian fit'
                        calcMethod = '2DGauss';
                    case 'max crossing profiles'
                        calcMethod = 'peakXTProfiles';
                    case 'estimate from event img'
                        calcMethod = 'estimate from event img';
                end
                % set show individual events to "on"
                showEventsFigOldVal = hObjs.check_showEventsFigs.Value;
                if showEventsFigOldVal < 1
                    hObjs.check_showEventsFigs.Value = 1;
                end
                % check if use normalized image or filtered raw
                if hObjs.check_useNormalizedImg.Value
                    img = imgData.imgDataXTfluoFN;
                else
                    % filter raw image
                    img = imgFiltering(imgData.imgDataXTfluoR, pxSzT, pxSzX);
                end
                findDetectedSparksParams(img, ...
                    detectedEvnts.detectedEvents(idx), ...
                    mainFig, calcMethod,...
                    detectedEvnts.detectedEventsRec(idx), ...
                    [], [], [], [], ...
                    hObjs.check_useNormalizedImg.Value); 
                % set back
                hObjs.check_showEventsFigs.Value = showEventsFigOldVal;           
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

    detectedEvnts.maskOfAcceptedSparks(idx) = state;
    detectedEvnts.typeOfEvent(idx) = typeOfEvent;

    % save changes
    setappdata(mainFig, 'sparkDetection', detectedEvnts)
    try
        calcSparkFreq(mainFig,true)
    catch
    end

catch 
   return
end    
    
end

