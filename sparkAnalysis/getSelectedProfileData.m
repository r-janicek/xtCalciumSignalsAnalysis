function prof = getSelectedProfileData(profNum,imgData,profileAnalysis,analysisType)

switch analysisType
    
    case 'spark recovery ryanodine'
        profData = profileAnalysis.selectedROIs.dataProfileRaw{profNum};
        [profData,~] = imgFiltering(profData,imgData.pxSzT,imgData.pxSzX);
        
    case 'spark recovery photolysis'
        %profData = profileAnalysis.selectedROIs.dataProfile{profNum};
        profData = profileAnalysis.selectedROIs.dataProfileRaw{profNum};

end

detectedEventsMask = profileAnalysis.selectedROIs.detectedEventsMask{profNum};
eventsPeaks = profileAnalysis.selectedROIs.eventsPeaks{profNum};
prevFitCoefSpRise = profileAnalysis.selectedROIs.prevFitCoef{profNum};

prof.t = imgData.t;
prof.y = mean(profData,1);

% set up baseline mask as percentile of all data 
% get fraction of data assigned to peaks 
peaksFr = sum(detectedEventsMask)/numel(detectedEventsMask);
baselineM = prof.y <= prctile( prof.y, round((1-peaksFr)*100) );
prof.baselineM = baselineM(:);
%prof.baselineM = ~detectedEventsMask;

prof.eventsM = detectedEventsMask;
prof.eventsPeaks = eventsPeaks;
prof.prevFitCoefSpRise = prevFitCoefSpRise;

end

