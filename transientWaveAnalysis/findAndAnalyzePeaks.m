function findAndAnalyzePeaks(~,~,analysisFig)

% get data
selectedEvent = getappdata(analysisFig,'selectedEvent');
hObjsA = getappdata(analysisFig,'hObjsA');
pxSzT = selectedEvent.pxSzT;
pxSzX = selectedEvent.pxSzX;
imgData = getappdata(analysisFig,'imgData');

% check if selected event was analyzed and accepted
try
    isAnalyzed = selectedEvent.analysis.accepted;
catch
    isAnalyzed = false;
end

switch selectedEvent.type
    
    case 'wave'
        y = selectedEvent.analysis.waveAnalysis.deskewedWaveProf;
        t = selectedEvent.ROIdata.dataROIs.t;
        x = selectedEvent.ROIdata.dataROIs.x;
        eventImg = selectedEvent.analysis.waveAnalysis.waveImg;
        eventImgDeskewed = selectedEvent.analysis.waveAnalysis.waveImgDeskewed;
        waveSpeed = selectedEvent.analysis.waveAnalysis.waveSpeed;
        tauD_regionWiseMed = median(cell2mat(selectedEvent.analysis.waveAnalysis.tauD_regionWise(:,1)));
        % axes to plot profile analysis
        hfBigAx = findall(0,'Tag','axesZoom');
        if hObjsA.check_expDPointSel.Value && ~isempty(hfBigAx)
            ax = findall(hfBigAx,'Type','axes');
        else
            ax = hObjsA.ax_deskewedProf;
        end
        
        axRes = hObjsA.ax_deskewedRes;
        if isAnalyzed
            splFit = selectedEvent.analysis.waveAnalysis.fitOfEvent.splFit;
            eventsParams = selectedEvent.analysis.waveAnalysis.fitOfEvent.calcParamsFromFit;
        end
        
    case 'transient'
        y = mean(selectedEvent.ROIdata.dataROIs.imgFN,1);
        t = selectedEvent.ROIdata.dataROIs.t;
        x = selectedEvent.ROIdata.dataROIs.x;
        waveSpeed = [];
        tauD_regionWiseMed = [];
        % axes to plot profile analysis
        % check zoom axes window exist
        hfBigAx = findall(0,'Tag','axesZoom');
        if hObjsA.check_expDPointSel.Value && ~isempty(hfBigAx)
            ax = findall(hfBigAx,'Type','axes');
        else
            ax = hObjsA.ax_orgImgProf;
        end
        
        axRes = hObjsA.ax_orgImgRes;
        if isAnalyzed
            splFit = selectedEvent.analysis.transientAnalysis.fitOfEvent.splFit;
            eventsParams = selectedEvent.analysis.transientAnalysis.fitOfEvent.calcParamsFromFit;
        end
        
    case 'caffeine' 
        y = selectedEvent.analysis.caffeineAnalysis.deskewedCaffeineProf;
        t = selectedEvent.ROIdata.dataROIs.t;
        x = selectedEvent.ROIdata.dataROIs.x;
        eventImg = selectedEvent.analysis.caffeineAnalysis.caffeineImg;
        eventImgDeskewed = selectedEvent.analysis.caffeineAnalysis.caffeineImgDeskewed;
        waveSpeed = [];
        tauD_regionWiseMed = median(cell2mat(selectedEvent.analysis.caffeineAnalysis.tauD_regionWise(:,1)));
        % axes to plot profile analysis
        % check zoom axes window exist
        hfBigAx = findall(0,'Tag','axesZoom');
        if hObjsA.check_expDPointSel.Value && ~isempty(hfBigAx)
            ax = findall(hfBigAx,'Type','axes');
        else
            ax = hObjsA.ax_deskewedProf;
        end
        axRes = hObjsA.ax_deskewedRes;
        if isAnalyzed
            splFit = selectedEvent.analysis.caffeineAnalysis.fitOfEvent.splFit;
            eventsParams = selectedEvent.analysis.caffeineAnalysis.fitOfEvent.calcParamsFromFit;
        end
        
end
  
if isAnalyzed
    if any(strcmp({'wave','caffeine'}, selectedEvent.type))
        % clear axes
        delete(hObjsA.ax_detEvent.Children)
        delete(hObjsA.ax_deskewed.Children)
        % show analysis
        image(eventImg, ...
            'YData',[min(x) max(x)],...
            'XData',[min(t) max(t)],...
            'CDataMapping','scaled', 'Parent',hObjsA.ax_detEvent);
        set(hObjsA.ax_detEvent, 'YGrid','off', 'FontSize',14)
        set(get(hObjsA.ax_detEvent,'Ylabel'), ...
            'String','x (\mum)', 'FontWeight','bold')
        set(get(hObjsA.ax_detEvent,'Xlabel'), ...
            'String','t (ms)', 'FontWeight','bold')
        set(get(hObjsA.ax_detEvent,'title'), ...
            'String','detected event', 'FontWeight','bold')
        
        image(eventImgDeskewed, ...
            'YData',[min(x) max(x)],...
            'XData',[min(t) max(t)],...
            'CDataMapping','scaled', 'Parent',hObjsA.ax_deskewed);
        set(hObjsA.ax_deskewed, 'XTick',[], 'YGrid','off', 'FontSize',14)
        set(get(hObjsA.ax_deskewed,'Ylabel'), ...
            'String','x (\mum)', 'FontWeight','bold')
        set(get(hObjsA.ax_deskewed,'title'), ...
            'String','deskewed event', 'FontWeight','bold')
        % set x axis
        hObjsA.ax_detEvent.XLim = [min(t) max(t)]; 
    end
    
    % delete previous fits and lines
    delete(findall(ax,'Tag','splineFit'))
    delete(findall(axRes,'Tag','splineFitRes'))
    delete(findall(ax,'Tag','fitLineRise'))
    delete(findall(ax,'Tag','profile'))
    
    % plot data and fit with spline
    line(t, y, 'Parent',ax, 'Color',[0 0 0 0.50], ...
        'LineStyle','-', 'LineWidth',1, 'Tag','profile');
    line(t, splFit, 'Parent',ax, 'Color',[1 0 0 0.75], ...
        'LineStyle','-', 'LineWidth',2, 'Tag','splineFit');
    line(t, splFit(:)-y(:), 'Parent',axRes, 'Color',[0 0 0], ...
        'LineStyle','-', 'LineWidth',1, 'Tag','splineFitRes');
else
    % two steps process, first fit with spline to smooth profile,
    % then find desired number of peaks
    % get spline fit params
    nK = str2double(hObjsA.h_edit_paramFit1.String);
    splOrd = str2double(hObjsA.h_edit_paramFit2.String);
    wholeEvntFit = fitOneWholeEvent(t, y, "spline", ...
        splineOrder = splOrd, ...
        spline_knots_seq = linspace(t(1), t(end), nK));
    splFit = wholeEvntFit.yFit;
    % find peaks
    N_events = str2double(hObjsA.h_edit_Npeaks.String);
    minPeakDist = str2double(hObjsA.h_edit_minPeakDist.String); % ms
    [pks, locs, w, p] = findpeaks(splFit, t, ...
        'SortStr','descend',...
        'MinPeakDistance',minPeakDist, ...
        'NPeaks',N_events);
    % sort peaks based on location
    [locs,I] = sort(locs);
    pks = pks(I);
    w = w(I);
    p = p(I);
    
    % delete previous fits and lines
    delete(findall(ax, 'Tag', 'fitLineRise'))
    
    % get selection of breaks positions for fitting 
    if isfield(selectedEvent,'expFitS') && N_events<2
        posOfSelPoints.expFitS = selectedEvent.expFitS;
        posOfSelPoints.expFitE = selectedEvent.expFitE;
        posOfSelPoints.bsFitS = selectedEvent.bsFitS;
        posOfSelPoints.bsFitE = selectedEvent.bsFitE;
        posOfSelPoints.peakFitPos = selectedEvent.peakFitPos;
        posOfSelPoints.peakFitVal = selectedEvent.peakFitVal;
        pks = selectedEvent.peakFitVal;
        locs = t(selectedEvent.peakFitPos);
    else
        posOfSelPoints = [];
    end
    
    % fit rise of event to get bs and t0
    bs_crit = str2double(hObjsA.h_edit_bsSens.String);
    smoothSpan = str2double(hObjsA.h_edit_knotsSpan.String);
    % parameters: [t0, tauR, A, bs]
    evntRiseFit = fitEventRise(t, y, "global", ...
        peaks_vals=pks, peaks_locs=locs, ax_prof=ax, ...
        smooth_span=smoothSpan, bs_crit=bs_crit, ...
        posOfSelPoints = posOfSelPoints, ...
        fittedLineColor=[0.2 1 0.2], ...
        fittedLineTag="fitLineRise");

    % delete previous fits and lines
    delete(findall(ax,'Tag','splineFit'))
    delete(findall(axRes,'Tag','splineFitRes'))
    delete(findall(ax,'Tag','profile'))
    
    % plot fit with spline
    line(t, y, 'Parent',ax, 'Color',[0 0 0 0.50], ...
        'LineStyle','-', 'LineWidth',1, 'Tag','profile');
    line(t, splFit, 'Parent',ax, 'Color',[1 0 0 0.75], ...
        'LineStyle','-', 'LineWidth',2, 'Tag','splineFit');
    line(t, splFit(:)-y(:), 'Parent',axRes, 'Color',[0 0 0], ...
        'LineStyle','-', 'LineWidth',1, 'Tag','splineFitRes');
    
    % put exp. rise fits on top
    uistack(evntRiseFit.h_line ,'up',2)
     
    % check if there is big axes window and also selection markers for expD
    % fitting
    if isfield(selectedEvent,'expFitS')
        % calculate phenomenological parameters of events from spline fit
        eventsParams = calcParametersOfGlobalEventsFromWholeProfileFit...
            (t, splFit, ...
            wholeEvntFit.x_ups, wholeEvntFit.yFit_ups, ...
            evntRiseFit.coef, evntRiseFit.startOfEvent, ...
            evntRiseFit.endOfEvent, N_events,...
            pks, locs, selectedEvent.expFitS, selectedEvent.expFitE);
        % move markers to top in axes
        uistack(findall(ax, 'Tag','expFitSelPoints'), 'top')     
    else
        % calculate phenomenological parameters of events from spline fit
        eventsParams = calcParametersOfGlobalEventsFromWholeProfileFit...
            (t, splFit, ...
            wholeEvntFit.x_ups, wholeEvntFit.yFit_ups, ...
            evntRiseFit.coef, evntRiseFit.startOfEvent, ...
            evntRiseFit.endOfEvent, N_events,...
            pks, locs, [], []);    
    end
    
end

% remove old
delete(findall(ax,'Tag','params'))
% show phenomenological parameters
try
    for i=2:size(eventsParams,1)
        
        line([eventsParams{i,1} eventsParams{i,4}], ...
            [eventsParams{i,3} eventsParams{i,3}],...
            'Parent',ax, 'Color','b', 'LineStyle','-', ...
            'LineWidth',2, 'Tag','params')
        line([eventsParams{i,5} eventsParams{i,6}], ...
            [eventsParams{i,7} eventsParams{i,7}],...
            'Parent',ax, 'Color','b', 'LineStyle','-', 'LineWidth',2,...
            'Tag','params')
        line([eventsParams{i,4} eventsParams{i,4}], ...
            [eventsParams{i,3} eventsParams{i,8}],...
            'Parent',ax, 'Color','b', 'LineStyle','-', 'LineWidth',2,...
            'Tag','params')
        line(eventsParams{i,18}(:,1), ...
            eventsParams{i,18}(:,2),...
            'Parent',ax, 'Color','b', 'LineStyle',':', 'LineWidth',2,...
            'Tag','params')  
    end
catch
end
% set up axes
set(ax, 'FontSize',14, 'YGrid','on')
set(get(ax,'Ylabel'), 'String','fluorescence (\DeltaF/F0)', ...
    'FontWeight','bold')
set(get(ax,'Xlabel'), 'String','t (ms)', 'FontWeight','bold')

set(axRes, 'YGrid','on', 'FontSize',14)
set(get(axRes,'Ylabel'), 'String','residuals', 'FontWeight','bold')
set(get(axRes,'Xlabel'), 'String','t (ms)', 'FontWeight','bold')

% set axes limits
ax.XLim = [min(t) max(t)];
ax.YLim = getAxisLimits([splFit(:);y(:)],5);
axRes.XLim = ax.XLim;

% parameters in text to show
if size(eventsParams,1)==2
    text(ax, t(end-5), max(y), ...
        {[sprintf('Amplitude: %0.2f',eventsParams{2,9}),'  \DeltaF/F_0'];...
        sprintf('TTP: %0.2f ms',eventsParams{2,10});...
        sprintf('FDHM: %0.2f ms',eventsParams{2,12});...
        sprintf('tauD fit: %0.2f ms',eventsParams{2,15});...
        sprintf('tauD median pxWise: %0.2f ms',tauD_regionWiseMed);...
        [sprintf('wave speed: %0.2f',waveSpeed),' \mum/s']},...
        'VerticalAlignment','top', 'Tag','params',...
        'HorizontalAlignment','right',...
        'FontSize',14)
end

if ~isAnalyzed
    % save data
    fitOfEvent.splFit = wholeEvntFit.yFit;
    fitOfEvent.calcParamsFromFit = eventsParams;
    
    switch selectedEvent.type
        case 'wave'
            selectedEvent.analysis.waveAnalysis.fitOfEvent = fitOfEvent; 
        case 'transient'
            selectedEvent.analysis.transientAnalysis.fitOfEvent = fitOfEvent;
        case 'caffeine'
            selectedEvent.analysis.caffeineAnalysis.fitOfEvent = fitOfEvent;   
    end
    setappdata(analysisFig,'selectedEvent',selectedEvent)
end

end

