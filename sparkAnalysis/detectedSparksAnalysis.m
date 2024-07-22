function detectedSparksAnalysis(h_O,~,mainFig)

% get data
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');
sparkDetection = getappdata(mainFig,'sparkDetection');

% selected method of parameters calculation
switch hObjs.calcSpParamMethodRBgroup.SelectedObject.String
    case '2D gaussian fit'
        calcMethod = '2DGauss';  
    case 'max crossing profiles'
        calcMethod = 'peakXTProfiles'; 
    case 'estimate from event img'
        calcMethod = 'estimate from event img';   
end

% first check if sparks were detected
if isempty(sparkDetection)
    if ~isfield(sparkDetection,'detectedEventsRec')
        errordlg('FIRST DETECT SPARKS!')
        return
    end
end

% remove old analysis
if isfield(sparkDetection,'eventParams')
    sparkDetection = rmfield(sparkDetection,'eventParams');
end
if isfield(sparkDetection,'maskOfAcceptedSparks')
    sparkDetection = rmfield(sparkDetection,'maskOfAcceptedSparks');
end

% check if use normalized image or filtered raw image to calculate spark
% parameters
if hObjs.check_useNormalizedImg.Value
    img = imgData.imgDataXTfluoFN;
else
    % filter raw image
    img = imgFiltering(imgData.imgDataXTfluoR, pxSzT, pxSzX);
end

% get sparks parameters
bsDetSensitivity = str2double(hObjs.h_bsDet_edit.String);
smoothSpan = str2double(hObjs.h_smooth_edit.String);
[sparkDetection.eventParams, sparkDetection.analyzedEvntsBrowserTbl] = ...
    findDetectedSparksParams( ...
    img, ...
    sparkDetection.detectedEvents, mainFig, ...
    calcMethod, sparkDetection.detectedEventsRec, ...
    [], [], [], [], hObjs.check_useNormalizedImg.Value, ...
    bsDetSensitivity, smoothSpan);

% save data
setappdata(mainFig, 'sparkDetection', sparkDetection);

% set numbers of events per page in case of saving them
if ~isempty(sparkDetection.analyzedEvntsBrowserTbl)
    possible_n_evnts_perFig = [1, 2, 4, 9, 12];
    [~, ind_n_evntsPerFig] = min( abs( possible_n_evnts_perFig - ...
        height(sparkDetection.analyzedEvntsBrowserTbl) ) );
    hObjs.h_edit_numOfEvntsPerPage.String = ...
        num2str(possible_n_evnts_perFig(ind_n_evntsPerFig));
end
% sparks filtering based on their parameters, also change colors of rectangles of rejected
% sparks
sparkParamsFiltering([], [], mainFig);

% show  individual calcium events
if hObjs.check_showEventsFigs.Value && ...
        ~isempty(sparkDetection.analyzedEvntsBrowserTbl)
    % close previous events viewer window, if any
    close(findall(0, 'Type','Figure', 'Tag','CaEventsBrowser'))
    % start events viewer window
    sparkDetection = getappdata(mainFig,'sparkDetection');
    eventsBrowserWindow(mainFig, sparkDetection.analyzedEvntsBrowserTbl)
end

end



