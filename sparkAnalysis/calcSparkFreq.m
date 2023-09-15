function calcSparkFreq(mainFig,corrected)
% calculate spark frequency in image

% data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

% calculate spark frequency
imgDataXTfluoFN = imgData.imgDataXTfluoFN;
pxSzX = imgData.pxSzX;
pxSzT = imgData.pxSzT;

sparkDetection = getappdata(mainFig,'sparkDetection');

% first check if sparks were detected
if isempty(sparkDetection)
    if ~isfield(sparkDetection,'detectedEventsRec')
        errordlg('FIRST DETECT SPARKS!')
        return
    end
end

% get parameters for spark frequency calculation
detectedEvents = sparkDetection.detectedEvents;
imgArea = (size(imgDataXTfluoFN,1)*pxSzX)*(size(imgDataXTfluoFN,2)*pxSzT/1000); % um*s

switch corrected
    case true
        maskOfAcceptedSparks = sparkDetection.maskOfAcceptedSparks;
        if isempty(maskOfAcceptedSparks)
            maskOfAcceptedSparks = true([numel(detectedEvents),1]);
        end
        sparkFreq = numel(detectedEvents(maskOfAcceptedSparks))*100/imgArea; % sparkFreq per 100um*s
        sparkDetection.correctedSparkFreq = sparkFreq;
        
        % set spark freq in text box in main window
        set(hObjs.txt_correctedSpF,'String',{sprintf('%0.1f',sparkFreq),['sp*100',char(956),'m-1*s-1' ]})
        
    otherwise
        sparkFreq = numel(detectedEvents)*100/imgArea; % sparkFreq per 100um*s
        sparkDetection.sparkFreq = sparkFreq;

        % set spark freq in text box in main window
        set(hObjs.txt_SpF,'String',{[sprintf('%0.1f  sp*100',sparkFreq),char(956),'m-1*s-1' ];'---'})
         
end

% save data
setappdata(mainFig,'sparkDetection',sparkDetection)

end

