function fitProfileBaseline(~, ~, fitFig)

set(fitFig,'Pointer','watch')
drawnow

hObjsFit = getappdata(fitFig,'hObjsFit');

% get selected profile data
prof = getappdata(fitFig,'selectedProf');

bsFit = fitBaseline(prof.t, prof.y, ...
    hObjsFit.popUpMenuBs.String{hObjsFit.popUpMenuBs.Value}, ...
    prof.baselineM, true, hObjsFit.ax_fit, ...
    [str2double(hObjsFit.h_edit_paramFitBs1.String), ...
     str2double(hObjsFit.h_edit_paramFitBs2.String)]);

% save baseline fit
prof.baselineFit = bsFit;
setappdata(fitFig,'selectedProf',prof);

% update baseline slider
bsSldVal = round(sum(prof.baselineM)/numel(prof.baselineM)*100);
if bsSldVal==100, bsSldVal=99; end
if bsSldVal==0, bsSldVal=1; end
hObjsFit.h_txt_sld_Bs.String = sprintf('mask of baseline as %d percentile',bsSldVal);
hObjsFit.sld_Bs.Value = bsSldVal;


set(fitFig,'Pointer','arrow')
drawnow

end

