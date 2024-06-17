function refitSelectedEvent(h_O, ~, ind, mainFig)
% refit selected event
% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
sparkDetection = getappdata(mainFig, 'sparkDetection');
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
% change values for baseline detection and smoothing
oldVal_bs = hObjs.h_bsDet_edit.String;
oldVal_smooth = hObjs.h_smooth_edit.String;
new_smooth_edit = getappdata(h_O.Parent.Parent, 'h_smooth_reFit');
new_bs_edit = getappdata(h_O.Parent.Parent, 'h_bsDet_reFit');
% set new values of baseline and smooth to main window
hObjs.h_bsDet_edit.String = new_bs_edit.String;
hObjs.h_smooth_edit.String = new_smooth_edit.String;

evntParams = findDetectedSparksParams(img, ...
    sparkDetection.detectedEvents(ind), ...
    mainFig, calcMethod,...
    sparkDetection.detectedEventsRec(ind), ...
    [], [], [], [], ...
    hObjs.check_useNormalizedImg.Value);
% keep the new figure, close the old one
close(h_O.Parent.Parent)
% update parameters of event in final structure
if isfield(sparkDetection, 'eventParams')
    params = fieldnames(sparkDetection.eventParams);
    for p = 1:numel(params)
        sparkDetection.eventParams.(params{p})(ind) = evntParams.(params{p});
    end
end
% save data
setappdata(mainFig, 'sparkDetection', sparkDetection);

% sparks filtering based on their parameters, also change colors of rectangles of rejected
% sparks
sparkParamsFiltering([], [], mainFig);

% change back smooth and baseline values in main window
hObjs.h_bsDet_edit.String = oldVal_bs;
hObjs.h_smooth_edit.String = oldVal_smooth;

end