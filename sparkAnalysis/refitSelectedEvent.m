function refitSelectedEvent(h_O, ~, mainFig, hf_evntBrowser)
% refit selected event
% get data
imgData = getappdata(mainFig,'imgData');
hObjs_evntsBrowser = getappdata(hf_evntBrowser, 'hObjs');
hObjs_mainFig = getappdata(mainFig, 'hObjs');
sparkDetection = getappdata(mainFig, 'sparkDetection');
analyzedEvntsBrowserTbl =  sparkDetection.analyzedEvntsBrowserTbl;
% selected event
ind = round(hObjs_evntsBrowser.sld_evnts.Value);
evnt = analyzedEvntsBrowserTbl(ind, :);

% check if use normalized image or filtered raw
if hObjs_mainFig.check_useNormalizedImg.Value
    img = imgData.imgDataXTfluoFN;
else
    % filter raw image
    img = imgFiltering(imgData.imgDataXTfluoR, imgData.pxSzT, imgData.pxSzX);
end
% refit event
[eventParams, refittedEvntTbl] = findDetectedSparksParams(img, ...
    sparkDetection.detectedEvents(ind), ...
    mainFig, evnt.calcMethod{1},...
    sparkDetection.detectedEventsRec(ind), ...
    [], [], [], [], ...
    hObjs_mainFig.check_useNormalizedImg.Value, ...
    str2double(hObjs_evntsBrowser.h_bsDet_reFit.String), ...
    str2double(hObjs_evntsBrowser.h_smooth_reFit.String));

% update parameters of event in final structure
if isfield(sparkDetection, 'eventParams')
    params = fieldnames(sparkDetection.eventParams);
    for p = 1:numel(params)
        sparkDetection.eventParams.(params{p})(ind) = eventParams.(params{p});
    end
end
% update table of analyzed events
refittedEvntTbl = [refittedEvntTbl, ...
    sparkDetection.analyzedEvntsBrowserTbl(ind, 'maskOfAcceptedSparks')];
sparkDetection.analyzedEvntsBrowserTbl(ind, :) = refittedEvntTbl;

% save data
setappdata(mainFig, 'sparkDetection', sparkDetection);

% sparks filtering based on their parameters, also change colors of rectangles of rejected
% sparks
sparkParamsFiltering([], [], mainFig);

% update events browser
showSelectedEventInBrowser(hObjs_evntsBrowser.sld_evnts, [], ...
    hf_evntBrowser, mainFig)

end