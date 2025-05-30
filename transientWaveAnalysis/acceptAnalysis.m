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
hf = figure('Name',sprintf('%s',evenROIname), 'units','normalized',...
    'Position',[0.1 0 0.4 1], 'Visible','off');
set(hf, 'PaperPositionMode', 'auto',...
    'PaperOrientation', 'portrait',...
    'PaperType', 'A4',...
    'Tag','imgsOfEventsFromAnalysis');

switch selectedEvent.type
    case 'transient'
        text_ax = axes('Parent',hf, ...
            'Units','normalized', 'Position',[0.1 0.94 0.8 0.03]);
        set(text_ax,'Visible','off')
        text((max(get(text_ax,'XLim'))-min(get(text_ax,'XLim')))/2,...
            min(get(text_ax,'YLim')), ...
            hObjsA.h_txt_fitNum.String,...
            'Parent',text_ax, 'FontUnits','normalized', ...
            'FontWeight','bold', 'FontSize',0.7, ...
            'VerticalAlignment','bottom', ...
            'HorizontalAlignment','center')
        % image and profile
        axImg = copyobj(hObjsA.ax_orgImg,hf);
        axProf = copyobj(hObjsA.ax_orgImgProf,hf);
        axImg.Position = [0.1 0.61 0.85 0.3];
        axProf.Position = [0.1 0.4 0.85 0.2];
        % add colorbar
        colormap(axImg,parula(256))
        h_img = findall(axImg,'Type','Image');
        clim(axImg, getAxisLimits(h_img.CData, 1))
        h_cb = colorbar(axImg,'west');
        h_cb.Position = [0.05 0.81 0.01 0.1];
        h_cb.Label.String = '\DeltaF/F_0';
        % show parameters of transients
        allTransParamsData = ...
            eventAnalysis.transientAnalysis.fitOfEvent.calcParamsFromFit(:,[9,10,11,12,14,15]);
        paramsNames = extractBefore(allTransParamsData(1,:)," (");
        paramsUnits = regexp(allTransParamsData(1,:), '[(].*[)]', 'match');
        paramsUnits = horzcat(paramsUnits{:});

        paramsNames = [{' '}, paramsNames];
        paramsUnits = [{'event #'}, paramsUnits];
        % replace amplitude to A and delta to Δ
        %paramsNames = strrep(paramsNames, 'amplitude', 'A');
        paramsUnits = strrep(paramsUnits, 'deltaF/F0', [char(916),'F/F0']);
        % replace _ with \_
        paramsNames = strrep(paramsNames, '_', '\_');
        paramsNums = [ num2cell((1:size(allTransParamsData(2:end,:),1)))', ...
            allTransParamsData(2:end,:) ];
        strTransientsParams = [ 
            paramsNames; ...
            paramsUnits; ...
            cellfun(@(x) num2str(round(x,2)), paramsNums, 'UniformOutput',0) ...
                ];        
        % pad strings
        for col = 1:size(strTransientsParams,2)
            n_charsToPad = max(8, ...
                max(cellfun(@numel,strTransientsParams(:,col),'UniformOutput',1)) );
            strTransientsParams(:,col) = pad(strTransientsParams(:,col), ...
                n_charsToPad, 'both', ' ');
        end
        strTransientsParams = [strTransientsParams(1:2,:); ...
                               repmat({'-'},size(paramsNames)); ...
                               strTransientsParams(3:end,:)];
        for col = 1:size(strTransientsParams,2)
            strTransientsParams(:,col) = pad(strTransientsParams(:,col), ...
                'both', '-');
        end
        % add empty space where \ was used
        for row = 1:size(strTransientsParams,1)
            for col = 1:size(strTransientsParams,2)
                n_ToAdd = count(strTransientsParams{row,col},'\');
                if n_ToAdd~=0
                    strTransientsParams{row,col} = ...
                        [strTransientsParams{row,col}, repmat(' ', [1 n_ToAdd])];
                end
            end
        end
        strTransientsParams = join(strTransientsParams, '|', 2);
        % show text
        text_ax2 = axes('Parent',hf, ...
            'Units','normalized', 'Position',[0.1 0.05 0.85 0.3], ...
            'Visible','off');
        text(min(get(text_ax2,'XLim')),...
            max(get(text_ax2,'YLim')), ...
            strTransientsParams,...
            'Parent',text_ax2, 'FontUnits','points', ...
            'FontWeight','normal', ...
            'FontSize',round(16*(get(0,'ScreenPixelsPerInch')/72)), ...
            'FontName','FixedWidth', ...
            'VerticalAlignment','top', ...
            'HorizontalAlignment','left')
              
    otherwise
        % invisible axes for text
        text_ax = axes('Parent',hf, ...
            'Units','normalized', 'Position',[0.1 0.94 0.8 0.03]);
        set(text_ax,'Visible','off')
        text((max(get(text_ax,'XLim'))-min(get(text_ax,'XLim')))/2,...
            min(get(text_ax,'YLim')), ...
            hObjsA.h_txt_fitNum.String,...
            'Parent',text_ax, 'FontUnits','normalized', ...
            'FontWeight','bold', 'FontSize',0.7, ...
            'VerticalAlignment','bottom', ...
            'HorizontalAlignment','center')

        axImg = copyobj(hObjsA.ax_orgImg,hf);
        axdetImg = copyobj(hObjsA.ax_detEvent,hf);
        axdeskwImg = copyobj(hObjsA.ax_deskewed,hf);
        axProf = copyobj(hObjsA.ax_deskewedProf,hf);
        axImg.Position = [0.1 0.71 0.85 0.2];
        axdetImg.Position = [0.1 0.49 0.85 0.2];
        axdetImg.XTick = [];
        axdetImg.XLabel = [];
        axdeskwImg.Position = [0.1 0.27 0.85 0.2];
        axProf.Position = [0.1 0.05 0.85 0.21];
        
        % add colorbar to first axes
        colormap(axImg,parula(256))
        h_img = findall(axImg,'Type','Image');
        clim(axImg,[ floor(prctile(h_img.CData(:),1)*10)/10 ...
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

