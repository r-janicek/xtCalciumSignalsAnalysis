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
sparkDetection.eventParams = findDetectedSparksParams( ...
    img, ...
    sparkDetection.detectedEvents, mainFig, ...
    calcMethod, sparkDetection.detectedEventsRec, ...
    [], [], [], [], hObjs.check_useNormalizedImg.Value);

% save data
setappdata(mainFig, 'sparkDetection', sparkDetection);

% sparks filtering based on their parameters, also change colors of rectangles of rejected
% sparks
sparkParamsFiltering([], [], mainFig);

end



