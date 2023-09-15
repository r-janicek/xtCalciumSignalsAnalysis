function fitSelectedEvent(~,~,analysisFig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

keyboard

% only fit data, 

set(analysisFig,'Pointer','watch')
drawnow

% get data
hObjsA = getappdata(analysisFig,'hObjsA');
% selected event number
currEventNum = hObjsA.sld_fitNum.Value;
% get current analysis
eventAnalysisAll = getappdata(analysisFig,'eventAnalysis');
eventAnalysis = eventAnalysisAll(currEventNum,:);

if ~isempty(regexp(eventAnalysis.roiName{1},'wave','once'))
    % wave
    keyboard
    t = eventAnalysis.waveAnalysis{1}.fitOfDeskewedWaveProf.t;
    y = eventAnalysis.waveAnalysis{1}.deskewedWaveProf;
    tauDpxWise = eventAnalysis.waveAnalysis{1}.tauD_pixelWise; 
    waveSpeed = eventAnalysis.waveAnalysis{1}.waveSpeed; 
    
    outFitEvent = fitEventProfile(t,y,analysisFig,hObjsA.ax_deskewedProf,hObjsA.ax_deskewedRes,tauDpxWise,waveSpeed);
    
    % save new fit
    eventAnalysis.waveAnalysis{1}.fitOfDeskewedWaveProf = outFitEvent;
    eventAnalysisAll(currEventNum,:) = eventAnalysis;

    
elseif ~isempty(regexp(eventAnalysis.roiName{1},'caffeine','once'))
    keyboard

    
else
    % transient
    t = eventAnalysis.transientAnalysis{1}.fitOfTransientsProf.t;
    y = eventAnalysis.transientAnalysis{1}.transientsProf;
    
    outFitEvent = fitEventProfile(t,y,analysisFig,hObjsA.ax_orgImgProf,hObjsA.ax_orgImgRes,[],[]);
    
    % save new fit
    eventAnalysis.transientAnalysis{1}.fitOfTransientsProf = outFitEvent;
    eventAnalysisAll(currEventNum,:) = eventAnalysis;
    
end

% save data  
setappdata(analysisFig,'eventAnalysis',eventAnalysisAll)

set(analysisFig,'Pointer','arrow')
drawnow

end

