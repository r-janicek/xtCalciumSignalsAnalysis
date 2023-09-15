function sparkParamsFiltering(~,~,mainFig)

% get data
hObjs = getappdata(mainFig,'hObjs');
sparkDetection = getappdata(mainFig,'sparkDetection');

% do filtering on sparks parameters, remove rectangles of refused one
if ~isfield(sparkDetection,'eventParams')
    return
end

% set detected sparks rectangles to only red color
arrayfun(@(x) set(x,'EdgeColor','r'), sparkDetection.detectedEventsRec )

% check amplitude
mA = sparkDetection.eventParams.amplitude > str2double(hObjs.h_edit_noise.String);

% check minimum FDHM
mFDHM = sparkDetection.eventParams.FDHM > str2double(hObjs.h_edit_spFDHM.String);

% check minimum FWHM
mFWHM = sparkDetection.eventParams.FWHM > str2double(hObjs.h_edit_spFWHM.String);

% filter detected sparks
maskOfAcceptedSparks = mA & mFDHM & mFWHM;
maskOfAcceptedSparks = maskOfAcceptedSparks';

% set accepted sparks rectangle color to red, not accepted to black
arrayfun(@(x) set(x,'EdgeColor','k'), sparkDetection.detectedEventsRec(~maskOfAcceptedSparks) )
arrayfun(@(x) set(x,'EdgeColor','k'), sparkDetection.detectedEventsMask(~maskOfAcceptedSparks) )

strRejectedSparks = arrayfun(@(x) get(x,'Tag'), ...
    sparkDetection.detectedEventsRec(~maskOfAcceptedSparks), ...
    'UniformOutput', 0 );
% get rectangle texts
rectTxt = findall(hObjs.ax_img,'Tag','detectedEventRecText');
rectTxtStr = arrayfun(@(x) get(x,'String'), rectTxt, 'UniformOutput', 0 );
% find mask for rextangle texts, it depends how findall is searching and
% sorting results

mStrRecTxt_rejectedSparks = ismember( str2num(char(rectTxtStr)), str2num(char(strRejectedSparks)) );

% set color of text, which belongs to rectangles
arrayfun(@(x) set(x,'Color','k'), rectTxt(mStrRecTxt_rejectedSparks) )

% save data
sparkDetection.maskOfAcceptedSparks = maskOfAcceptedSparks;
setappdata(mainFig,'sparkDetection',sparkDetection)

% calculate corrected spark frequency
calcSparkFreq(mainFig,true)

end

