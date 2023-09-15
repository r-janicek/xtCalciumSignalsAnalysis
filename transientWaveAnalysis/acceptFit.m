function acceptFit(~,~,analysisFig)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% get data
hObjsA = getappdata(analysisFig,'hObjsA');
% selected event number
currEventNum = hObjsA.sld_fitNum.Value;
% get current analysis
eventAnalysisAll = getappdata(analysisFig,'eventAnalysis');
eventAnalysis = eventAnalysisAll(currEventNum,:);
eventAnalysis.acceptedFitOfEvent = true;
eventAnalysisAll(currEventNum,:) = eventAnalysis;

% save data to mainFig
mainFig = getappdata(analysisFig,'mainFig');
selectedROIs = getappdata(mainFig,'selectedROIs');
selectedROIs.analysis(strcmp(selectedROIs.roiName,eventAnalysis.roiName{1})) = {eventAnalysis};

% save data
setappdata(analysisFig,'eventAnalysis',eventAnalysisAll)
setappdata(mainFig,'selectedROIs',selectedROIs)

% set up window
% disable baseline fitting
hObjsA.popUpMenuEventsFcn.Enable = 'off';
hObjsA.h_edit_paramFit1.Enable = 'off';
hObjsA.h_edit_paramFit2.Enable = 'off';
hObjsA.h_edit_Npeaks.Enable = 'off';
hObjsA.check_pieceWise.Enable = 'off';
hObjsA.check_signOfPeak.Enable = 'off';

% disable events fitting
hObjsA.h_pb_fit.Enable = 'off';
hObjsA.h_pb_acceptFit.Enable = 'off';
hObjsA.check_fitWeights.Enable = 'off';
hObjsA.h_pb_stopFiting.Enable = 'off';

end

