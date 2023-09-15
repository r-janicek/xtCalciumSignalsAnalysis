function acceptFit(~,~,fitFig)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% get data
hObjsFit = getappdata(fitFig,'hObjsFit');
prof = getappdata(fitFig,'selectedProf');
profileAnalysis = getappdata(fitFig,'profileAnalysis');
selectedROIs = profileAnalysis.selectedROIs;

% get profile number
profnum = hObjsFit.sld_fitNum.Value;

% height of selected rois table
hT = height(selectedROIs);

if any(strcmp(selectedROIs.Properties.VariableNames,'wholeProfileFit'))
   
    selectedROIs.wholeProfileFit{profnum} = prof;
  
else
    
    wholeProfileFit = num2cell(nan(hT,1));
    wholeProfileFit(profnum) = {prof};
        
    selectedROIs.wholeProfileFit = wholeProfileFit; 
              
end

% save data to main fig
profileAnalysis.selectedROIs = selectedROIs;
setappdata(fitFig,'profileAnalysis',profileAnalysis)

mainFig = getappdata(fitFig,'mainFig');

setappdata(mainFig,'profileAnalysis',profileAnalysis)

% set up window
% disable baseline fitting
hObjsFit.popUpMenuBs.Enable = 'off';
hObjsFit.h_edit_paramFitBs1.Enable = 'off';
hObjsFit.h_edit_paramFitBs2.Enable = 'off';
hObjsFit.h_pb_fitBs.Enable = 'off';
hObjsFit.h_pb_normProf.Enable = 'off';
hObjsFit.h_pb_undoNormProf.Enable = 'off';

hObjsFit.maskButtonGroup.SelectedObject = hObjsFit.rbutton2;
hObjsFit.rbutton1.Enable = 'off';
hObjsFit.rbutton2.Enable = 'on';

% disable events fitting
hObjsFit.h_pb_fit.Enable = 'off';
hObjsFit.h_pb_acceptFit.Enable = 'off';

end

