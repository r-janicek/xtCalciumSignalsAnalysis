function fitNum_slider(h_sld,~,analysisFig)

%% get data
mainFig = getappdata(analysisFig,'mainFig');
selectedROIs = getappdata(mainFig,'selectedROIs');
hObjsA = getappdata(analysisFig,'hObjsA');
imgData = getappdata(analysisFig,'imgData');

% selected event number
currEventNum = hObjsA.sld_fitNum.Value;

% create selectedEvent structure
selectedEvent.ROIdata = selectedROIs(currEventNum,:);
selectedEvent.pxSzT = imgData.pxSzT;
selectedEvent.pxSzX = imgData.pxSzX;
roiName = selectedEvent.ROIdata.roiName{1};
% save type of event
if ~isempty(regexp(roiName,'wave','once'))
    selectedEvent.type = 'wave';
elseif ~isempty(regexp(roiName,'caffeine','once'))
    selectedEvent.type = 'caffeine';
else
    selectedEvent.type = 'transient';
end
% check if selected event was analyzed and accepted
try
    isAnalyzed = selectedEvent.ROIdata.analysis{1}.accepted;
catch
    isAnalyzed = false;
end
% save analysis
if isAnalyzed
    selectedEvent.analysis = selectedEvent.ROIdata.analysis{1};
end

% save selected event
setappdata(analysisFig,'selectedEvent',selectedEvent)

% set up window for fitting
hObjsA.h_edit_paramFit1.String = ...
    num2str( ceil(max(selectedEvent.ROIdata.dataROIs.t)/50) );

% set number of peaksfor a new event to 1
hObjsA.h_edit_Npeaks.String = num2str(1); 


%% show filtered and normalized image of selected event
% clear axes
arrayfun(@(x) cla(x),findobj(analysisFig,'Type','axes'))

% show image
image(selectedEvent.ROIdata.dataROIs.imgFN,...
    'YData',[min(selectedEvent.ROIdata.dataROIs.x) max(selectedEvent.ROIdata.dataROIs.x)],...
    'XData',[min(selectedEvent.ROIdata.dataROIs.t) max(selectedEvent.ROIdata.dataROIs.t)],...
    'CDataMapping','scaled','Parent',hObjsA.ax_orgImg);
set(hObjsA.ax_orgImg,'XTick',[],'YGrid','off','FontSize',14)
set(get(hObjsA.ax_orgImg,'Ylabel'),'String','x (\mum)','FontWeight','bold')
set(get(hObjsA.ax_orgImg,'title'),'String','image filtered and normalized','FontWeight','bold')


%% analyze image or show already analyzed image and accepted analysis
if isAnalyzed
    % just show resutls, use analyze peaks function without fitting
    findAndAnalyzePeaks([],[],analysisFig)
    
else 
    % do whole analysis
    analyzeEvent([],[],analysisFig);
    
end
   

%% set up analysis window
hObjsA.h_txt_fitNum.String = sprintf('analysis of event #%d:  %s',h_sld.Value,roiName);
hObjsA.h_txt_fitNumSld.String = sprintf('event #%d;  %s',h_sld.Value,roiName);

if isAnalyzed
    % set up window
    % disable fitting
    hObjsA.h_edit_paramFit1.Enable = 'off';
    hObjsA.h_edit_paramFit2.Enable = 'off';
    hObjsA.h_edit_Npeaks.Enable = 'off';
    hObjsA.check_signOfPeak.Enable = 'off';
    hObjsA.h_edit_bsSens.Enable = 'off';
    
    % disable mask correction
    hObjsA.h_edit_ellipseAxes1.Enable = 'off';
    hObjsA.h_edit_ellipseAxes2.Enable = 'off';
    hObjsA.h_edit_treshDet.Enable = 'off';
    hObjsA.h_edit_filtSzT.Enable = 'off';
    hObjsA.h_edit_filtSzX.Enable = 'off';
    hObjsA.h_edit_knotsSpan.Enable = 'off';
    hObjsA.h_pb_createROIcircle.Enable = 'off';
    
    % disable pushbuttons
    hObjsA.h_pb_acceptAnalysis.Enable = 'off';
    hObjsA.h_pb_refit.Enable = 'off';
    hObjsA.h_pb_detect.Enable = 'off';
    hObjsA.h_pb_deskew.Enable = 'off';
    hObjsA.h_pb_eventAnalysis.Enable = 'off';
    
else
    % set up window
    % enable fitting
    hObjsA.h_edit_paramFit1.Enable = 'on';
    hObjsA.h_edit_paramFit2.Enable = 'on';
    hObjsA.h_edit_Npeaks.Enable = 'on';
    hObjsA.check_signOfPeak.Enable = 'on';
    hObjsA.h_edit_bsSens.Enable = 'on';
    
    % enable mask correction
    hObjsA.h_edit_ellipseAxes1.Enable = 'on';
    hObjsA.h_edit_ellipseAxes2.Enable = 'on';
    hObjsA.h_edit_treshDet.Enable = 'on';
    hObjsA.h_edit_filtSzT.Enable = 'on';
    hObjsA.h_edit_filtSzX.Enable = 'on';
    hObjsA.h_edit_knotsSpan.Enable = 'on';
    hObjsA.h_pb_createROIcircle.Enable = 'on';
    
    % enable pushbuttons
    hObjsA.h_pb_acceptAnalysis.Enable = 'on';
    hObjsA.h_pb_refit.Enable = 'on';
    hObjsA.h_pb_detect.Enable = 'on';
    hObjsA.h_pb_deskew.Enable = 'on';
    hObjsA.h_pb_eventAnalysis.Enable = 'on';

end

end
