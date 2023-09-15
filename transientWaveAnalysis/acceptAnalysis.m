function acceptAnalysis(~,~,analysisFig)

% get data
hObjsA = getappdata(analysisFig,'hObjsA');

% selected event number
currEventNum = hObjsA.sld_fitNum.Value;

% get current analysis
selectedEvent = getappdata(analysisFig,'selectedEvent');
eventAnalysis = selectedEvent.analysis;
eventAnalysis.accepted = true;
eventAnalysis.type = selectedEvent.type;
evenROIname = selectedEvent.ROIdata.roiName{1};

% save data to mainFig
mainFig = getappdata(analysisFig,'mainFig');
selectedROIs = getappdata(mainFig,'selectedROIs');
selectedROIs.analysis(strcmp(selectedROIs.roiName,evenROIname)) = {eventAnalysis};

% save data
setappdata(mainFig,'selectedROIs',selectedROIs)

%% copy images of individual analysis to invisible figure
hf = figure('Name',sprintf('%s',evenROIname),'units','normalized',...
    'Position',[0.1 0 0.4 1],'Visible','off');
set(hf, 'PaperPositionMode', 'auto',...
    'PaperOrientation', 'portrait',...
    'PaperType', 'A4',...
    'Tag','imgsOfEventsFromAnalysis');

%uistack(hf,'bottom')
%hf.Visible = 'off';

switch selectedEvent.type
    case 'transient'
        titleStr = copyobj(hObjsA.h_txt_fitNum,hf);
        axImg = copyobj(hObjsA.ax_orgImg,hf);
        axProf = copyobj(hObjsA.ax_orgImgProf,hf);
        
        titleStr.Position = [0.1 0.94 0.8 0.03];
        axImg.Position = [0.1 0.61 0.85 0.3];
        axProf.Position = [0.1 0.4 0.85 0.2];
        % add colorbar
        colormap(axImg,parula(256))
        h_img = findall(axImg,'Type','Image');
        caxis(axImg,[ floor(prctile(h_img.CData(:),1)*10)/10 ...
            ceil(prctile(h_img.CData(:),99.9)) ])
        h_cb = colorbar(axImg,'west');
        h_cb.Position = [0.05 0.81 0.01 0.1];
        h_cb.Label.String = '\DeltaF/F_0';
        
        uit = uitable(hf,...
            'Data',eventAnalysis.transientAnalysis.fitOfEvent.calcParamsFromFit(2:end,[1,3,8,9,10,11,12]),...
            'units','normalized',...
            'Position',[0.1 0.05 0.85 0.3]);
        uit.ColumnName = cellfun(@(x) sprintf('<HTML> <font size="4"> <b> %s </b> </HTML>',x),...
            eventAnalysis.transientAnalysis.fitOfEvent.calcParamsFromFit(1,[1,3,8,9,10,11,12]),...
            'UniformOutput',0);
        uit.ColumnEditable = false(1,7);
        uit.ColumnFormat = {'short g','short g','short g','short g',...
            'short g','short g','short g'};
        uit.FontSize = 14;
            
    otherwise
        titleStr = copyobj(hObjsA.h_txt_fitNum,hf);
        axImg = copyobj(hObjsA.ax_orgImg,hf);
        axdetImg = copyobj(hObjsA.ax_detEvent,hf);
        axdeskwImg = copyobj(hObjsA.ax_deskewed,hf);
        axProf = copyobj(hObjsA.ax_deskewedProf,hf);
        
        titleStr.Position = [0.1 0.94 0.8 0.03];
        axImg.Position = [0.1 0.71 0.85 0.2];
        axdetImg.Position = [0.1 0.49 0.85 0.2];
        axdetImg.XTick = [];
        axdetImg.XLabel = [];
        axdeskwImg.Position = [0.1 0.27 0.85 0.2];
        axProf.Position = [0.1 0.05 0.85 0.21];
        
        % add colorbar to first axes
        colormap(axImg,parula(256))
        h_img = findall(axImg,'Type','Image');
        caxis(axImg,[ floor(prctile(h_img.CData(:),1)*10)/10 ...
            ceil(prctile(h_img.CData(:),99.9)) ])
        h_cb = colorbar(axImg,'west');
        h_cb.Position = [0.05 0.81 0.01 0.1];
        h_cb.Label.String = '\DeltaF/F_0';
        
               
end


%% set up window
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

end

