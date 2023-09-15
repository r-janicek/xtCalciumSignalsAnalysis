function [eventParams, pairedEventsFit] = findDetectedEventsParamsPair_2DGaussFit(...
    img,statEvents,startOfEvent,endOfEvent,prevFitCoef,...
    posOfROI,mainFig,selectedROI)

% img = not normalized data, raw data
% statEvents = statistic from event detection
% startOfSpark,endOfSpark,prevFitCoef = from previous fitting in spark recovery analysis 
% posOfROI = from spark recovery analysis
if isempty(posOfROI)
    posOfROI = nan;
end

% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
% get results from previous fitting
wholeProfileFitResults = selectedROI.wholeProfileFit{1}.profFit.coefficientsOfFittedEvents;
wholeProfileFitResults = cell2mat(wholeProfileFitResults(2:end,2:end));
bsWholeProfFit = mean(selectedROI.wholeProfileFit{1}.baselineFit);
% mask of events
maskEvents = selectedROI.wholeProfileFit{1}.eventsM;

pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;
scrRes = get(0,'ScreenSize');

% % % get normalized image of selected area of profile
% % whImgR = imgData.imgDataXTfluoR;
% % whImgFN = imgData.imgDataXTfluoFN;
% %
% % % find raw img data in whole image, use cross-correlation
% % nimg = whImgR-mean(mean(whImgR));
% % imgR = img;
% % crr = xcorr2(nimg,imgR);
% % [ssr,snd] = max(crr(:));
% % [ij,ji] = ind2sub(size(crr),snd);
% % 
% % img = whImgFN(ij:-1:ij-size(img,1)+1,ji:-1:ji-size(img,2)+1);
% % img = rot90(img,2);

% find position of TPP in wholeProfile area
% in pixels
r_dataWholeArea = (size(selectedROI.dataWholeArea{1},1)-1)/2;
midPoint = selectedROI.positionOfRoi;
posOfTTP_inX_dataWholeArea = round( imgData.posOfTPPinScanLine - ...
    (imgData.cropROIpos(2) + midPoint-r_dataWholeArea - 1));

% raw data image
imgR = img;

% filter raw image
img = imgFiltering(img,pxSzT,pxSzX);

h_d = str2double(get(hObjs.h_edit_averageWidth,'String'));

n_px_t = ceil(h_d/pxSzX);
if mod(n_px_t,2)==0
    n_px_t = n_px_t - 1;
end

n_px_x = ceil(5/pxSzT);
if mod(n_px_x,2)==0
    n_px_x = n_px_x - 1;
end

% add waitbar
hw = waitbar(0,'0 %','Name','analyzing detected sparks...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(hw,'canceling',0)

N_pairs = length(statEvents)/2;

% estimate parameters from pairs, process whole pair 
for i = 1:N_pairs
    
    % model used to fit data
    switch selectedROI.wholeProfileFit{1}.profFit.eventModel
        case 'CaSpikeFun'
            gauss2Dmodel = 'tCaSpike_xEMG';
        otherwise  
            gauss2Dmodel = '2D_EMGxt';  % '2D_EMGxt'
    end     

    % first do single events in pair to get initial parameters for pair fitting  
    separateEventsFits = findDetectedEventsParamsLocalFun(...
        imgR-mean(imgR(:)),statEvents(2*i-1:2*i),...
        startOfEvent(2*i-1:2*i),...
        endOfEvent(2*i-1:2*i),...
        prevFitCoef((2*i-1:2*i),:),...
        mainFig,'2DGauss',gauss2Dmodel);
   
    % build mask of paired event in whole image
    imgMask = false(size(img));
    % add 1. event 
    imgMask(statEvents(2*i-1).SubarrayIdx{1},statEvents(2*i-1).SubarrayIdx{2}) = ...
        statEvents(2*i-1).Image;
    imgMaskFirst = imgMask;
    % add 2. event
    imgMask(statEvents(2*i).SubarrayIdx{1},statEvents(2*i).SubarrayIdx{2}) = ...
        statEvents(2*i).Image;
    imgMaskSecond = false(size(img));
    imgMaskSecond(statEvents(2*i).SubarrayIdx{1},statEvents(2*i).SubarrayIdx{2}) = ...
        statEvents(2*i).Image; 
     
    % get data for pair of events
    rowsOfEventsInPair = [statEvents(2*i-1).SubarrayIdx{1},statEvents(2*i).SubarrayIdx{1}];
    colsOfEventsInPair = [statEvents(2*i-1).SubarrayIdx{2},statEvents(2*i).SubarrayIdx{2}];
    
    pairRows = (min(rowsOfEventsInPair):max(rowsOfEventsInPair));
    pairCols = (min(colsOfEventsInPair):max(colsOfEventsInPair));
    
    % get starting and ending points of pair from whole profile fitting 
    weightedCentroidEvent_col_1 = round(statEvents(2*i-1).WeightedCentroid(1));
    weightedCentroidEvent_col_2 = round(statEvents(2*i).WeightedCentroid(1));
    %midPointPair = round(median(pairCols));
    
    % find start of event form ask of events from whole profile fitting
    sPairWholeProfFit = weightedCentroidEvent_col_1 - ...
        find(flipud(maskEvents(1:weightedCentroidEvent_col_1))==0,1,'first') + 2;
    
    if isempty(sPairWholeProfFit), sPairWholeProfFit=1; end
    ePairWholeProfFit = weightedCentroidEvent_col_2 + ...
        find(maskEvents(weightedCentroidEvent_col_2:end)==0,1,'first') - 2;
    if isempty(ePairWholeProfFit), ePairWholeProfFit=pairCols(end); end
    if ePairWholeProfFit<=sPairWholeProfFit, ePairWholeProfFit=pairCols(end); end
    try
        % check if positions of all two PP pulses are covered in image of the pair
        % check for closest
        if imgData.s_TPP(i*2-1) > sPairWholeProfFit
            pairCols = sPairWholeProfFit:ePairWholeProfFit;                  
        else
            % check what is wrong
            % try to find closest relevant TPP
            d_TPP = imgData.s_TPP-sPairWholeProfFit;
            [~,idx_d_TPP] = min(abs(d_TPP));
            if d_TPP(idx_d_TPP) > 0
                % most likely some triggered events were not taken into
                % analysis
                pairCols = sPairWholeProfFit:ePairWholeProfFit;
            else
                % 2PP is not included, do not use mask from whole profile
                % fit
                if ~isempty(startOfEvent)
                    sPair = startOfEvent(2*i-1);
                    ePair = endOfEvent(2*i);
                    if sPair>pairCols(1), sPair = pairCols(1); end
                    pairCols = sPair:ePair;
                end
               
            end 
            
        end
    catch
        if ~isempty(startOfEvent)
            sPair = startOfEvent(2*i-1);
            ePair = endOfEvent(2*i);
            if sPair>pairCols(1), sPair = pairCols(1); end
            pairCols = sPair:ePair;
        end
    end

   
    % adjusted event data, pair 
    imgPair = img(pairRows,pairCols);
    % mask of pair of events
    imgPair_m = imgMask(pairRows,pairCols);
    % masks of individual events
    imgPair_events_m = false(numel(pairRows),numel(pairCols));
    imgPair_events_m(:,:,1) = imgMaskFirst(pairRows,pairCols);
    imgPair_events_m(:,:,2) = imgMaskSecond(pairRows,pairCols);
    
    % create time and spatial axes
    t = (1:1:size(imgPair,2)).*pxSzT - pxSzT; % x axes = time
    x = (1:1:size(imgPair,1)).*pxSzX - pxSzX; % y axes = spatial
    t = t(:);
    x = x(:);
   
    % upscale 10x
    uspFactor = 10;
    t_ups = linspace(t(1),t(end),numel(t)*uspFactor); % x axes = time
    x_ups = linspace(x(1),x(end),numel(x)*uspFactor); % y axes = spatial
    t_ups = t_ups(:);
    x_ups = x_ups(:);
    
    % find row where is maximum
    % smoothed spark data with spline
    Ds = csaps({linspace(1,size(imgPair,1),size(imgPair,1)),linspace(1,size(imgPair,2),size(imgPair,2))},...
        imgPair,0.25,...
        {linspace(1,size(imgPair,1),size(imgPair,1)),linspace(1,size(imgPair,2),size(imgPair,2))});
    
    % mask of pixels higher than 90. percentile
    m_spD = Ds>prctile(Ds(:),90);
    CC_event = bwconncomp(m_spD,8);
    
    %calculate properties of all CC (regions) in spark area, to get centre
    %of spark
    statOfSubRegions = regionprops(CC_event,Ds, 'WeightedCentroid','Area','SubarrayIdx');
    
    % find biggest region
    [~,p] = max([statOfSubRegions.Area]);
    
    % get position of centre of spark
    %centr = statEvents(i).WeightedCentroid;
    %r_m = round(centr(2)) - rows_e(1)+1;
    center = statOfSubRegions(p).WeightedCentroid;
    r_m = round(center(2));
    
    % calculate parameters, fit spark data with 2D exponentially modified gaussian
    % create time and spatial grids
    [T,X] = meshgrid(t,x);
    [T_ups,X_ups] = meshgrid(t_ups,x_ups);
    
    % get parameters of individual fits of event in pair from previous whole profile fitting
    event1_fitParams = wholeProfileFitResults(2*i-1,:);
    event2_fitParams = wholeProfileFitResults(2*i,:);
     
    % initial parameters from separate events fitting
    switch gauss2Dmodel
        case 'tCaSpike_xEMG'
            % 2D CaSpike function in time and expModGauss in spatial direction
            %coefficients: [t0, FM, tA, tI, FI, m_x, sd_x, expTauD_x]
  
            % estimate of maxima of events from individual 2D fits
            maxOfEvents = [max(separateEventsFits{1}.dataEventFit(:)), ...
                max(separateEventsFits{2}.dataEventFit(:))];
  
            % get estimated position of maximum of events, take weighted
            % centroid position
            % adjust for beggining of ROI, time starts at 0 
            posOfEvents_t = ([statEvents(2*i-1).WeightedCentroid(1),...
                statEvents(2*i).WeightedCentroid(1)] - pairCols(1) + 1).*pxSzT;
            
            % from previous whole profile fitting [t0, FM, tA, tI, FI]
            % get estimate of t0 for events from previous fitting
            t0_prevFit = [event1_fitParams(1),event2_fitParams(1)] - pairCols(1)*pxSzT;
            for ii=1:numel(t0_prevFit)
                [~,ind_t0] = min( abs(t-t0_prevFit(ii)) );
                t0_est(ii) = t(ind_t0);
            end
            
            % in um
            posOfEvents_x = [ separateEventsFits{1}.coef2DGauss(6) + ...
                abs(min(pairRows)-min(statEvents(2*i-1).SubarrayIdx{1}))*pxSzX, ...
                separateEventsFits{2}.coef2DGauss(6) + ...
                abs(min(pairRows)-min(statEvents(2*i).SubarrayIdx{1}))*pxSzX];
            sd_x = [separateEventsFits{1}.coef2DGauss(7), ...
                separateEventsFits{2}.coef2DGauss(7)];
            expTauD_x = [separateEventsFits{1}.coef2DGauss(8), ...
                separateEventsFits{2}.coef2DGauss(8)];
        
            prevFitParams.maxOfEvents = maxOfEvents;
            prevFitParams.posOfEvents_t = posOfEvents_t;
            prevFitParams.t0_est = t0_est;
            prevFitParams.posOfEvents_x = posOfEvents_x;
            prevFitParams.sd_x = sd_x;
            prevFitParams.expTauD_x = expTauD_x;
            prevFitParams.FM = [separateEventsFits{1}.coef2DGauss(2),...
                separateEventsFits{2}.coef2DGauss(2)];
            prevFitParams.tA = [separateEventsFits{1}.coef2DGauss(3),...
                separateEventsFits{2}.coef2DGauss(3)];
            prevFitParams.tI = [separateEventsFits{1}.coef2DGauss(4),...
                separateEventsFits{2}.coef2DGauss(4)];
            prevFitParams.FI = [separateEventsFits{1}.coef2DGauss(5),...
                separateEventsFits{2}.coef2DGauss(5)];
            
            
        case '2D_EMGxt'
            % 2D EMGtx
            % coefficients: [A, peakPos_t, sd_t, expTauD_t, peakPos_x, sd_x, expTauD_x, F0]
            maxOfEvents = [separateEventsFits{1}.coef2DGauss(1), ...
                separateEventsFits{2}.coef2DGauss(1)];
            
%             maxOfEvents = [max(separateEventsFits{1}.dataEventFit(:)), ...
%                 max(separateEventsFits{2}.dataEventFit(:))];
            
            % from previous whole profile fitting {'A','m','sd','tauD'}
            % in ms
             % adjust for beggining of ROI, time starts at 0 
            posOfEvents_t = [event1_fitParams(2),event2_fitParams(2)]- pairCols(1)*pxSzT;
            sd_t = [separateEventsFits{1}.coef2DGauss(3), ...
                separateEventsFits{1}.coef2DGauss(3)];
            expTauD_t = [separateEventsFits{1}.coef2DGauss(4), ...
                separateEventsFits{1}.coef2DGauss(4)];
            
            % in um
            posOfEvents_x = [ separateEventsFits{1}.coef2DGauss(5) + ...
                abs(min(pairRows)-min(statEvents(2*i-1).SubarrayIdx{1}))*pxSzX, ...
                separateEventsFits{2}.coef2DGauss(5) + ...
                abs(min(pairRows)-min(statEvents(2*i).SubarrayIdx{1}))*pxSzX];
            sd_x = [separateEventsFits{1}.coef2DGauss(6), ...
                separateEventsFits{2}.coef2DGauss(6)];
            expTauD_x = [separateEventsFits{1}.coef2DGauss(7), ...
                separateEventsFits{2}.coef2DGauss(7)];
            
            prevFitParams.maxOfEvents = maxOfEvents;
            prevFitParams.posOfEvents_t = posOfEvents_t;
            prevFitParams.sd_t = sd_t;
            prevFitParams.expTauD_t = expTauD_t;
            prevFitParams.posOfEvents_x = posOfEvents_x;
            prevFitParams.sd_x = sd_x;
            prevFitParams.expTauD_x = expTauD_x;
    end

    % get lower bounds for position of second event in pair
    % from 2PP pulses or position of first event
    if isfield(imgData,'s_TPP') && ~isempty(imgData.s_TPP)
        
        sPP_all = imgData.s_TPP;
        ePP_all = imgData.e_TPP;
        sPP = sPP_all(ismember(sPP_all,pairCols));
        ePP = ePP_all(ismember(ePP_all,pairCols));
        
        sPP_pair = sPP.*pxSzT;
        ePP_pair = ePP.*pxSzT;
        
    else
        % use positions of previous peaks as lower bound
        sPP_pair = [0,prevFitParams.posOfEvents_t(1:end-1)];
        %sPP_n = zeros(size(pos_peak));
    end
    
    % adjust for beggining of ROI, time starts at 0 
    sPP_pair = sPP_pair - pairCols(1)*pxSzT - pxSzT;
    ePP_pair = ePP_pair - pairCols(1)*pxSzT - pxSzT;
    % position of TPP in x axis in pairEvents area
    posOfTTP_inPairEventsArea = posOfTTP_inX_dataWholeArea - pairRows(1) - 1;
    posOfTTP_inPairEventsArea = posOfTTP_inPairEventsArea*pxSzX; % in um
 
    % try to fit multiple 2D ExpModGauss
    try       
        fit2DGauss = multiple2DGaussianFit( ...
            T,X,T_ups,X_ups,imgPair,imgPair_m,imgPair_events_m,...
            prevFitParams,sPP_pair,...
            gauss2Dmodel,bsWholeProfFit);
        
        D_fit_events = fit2DGauss.individualEventsFits;
        D_fit_ups_events = fit2DGauss.ups.individualEventsFits;

        % get parameters of individual events in pair
        t_fit_prof = zeros(size(D_fit_ups_events{2,1},2),numel(maxOfEvents));
        x_fit_prof = zeros(size(D_fit_ups_events{2,1},1),numel(maxOfEvents));
        for p = 1:numel(maxOfEvents)
            
            % get t and x profiles crossing max of fit
            [rFitMax(p),cFitMax(p)] = find(D_fit_events{2,p} == max(D_fit_events{2,p}(:)),1);
            % t_data_prof = D(rFitMax,:); t_data_prof = t_data_prof(:);
            % x_data_prof = D(:,cFitMax); x_data_prof = x_data_prof(:);
            [rFitMaxUps(p),cFitMaxUps(p)] = ...
                find(D_fit_ups_events{2,p} == max(D_fit_ups_events{2,p}(:)),1);
            t_fit_prof(:,p) = D_fit_ups_events{2,p}(rFitMaxUps(p),:)+fit2DGauss.baselineEst;
            x_fit_prof(:,p) = D_fit_ups_events{2,p}(:,cFitMaxUps(p))+fit2DGauss.baselineEst;

            eP_2Dfit(p,1) = getParametersOfEventProfile(t_ups,t_fit_prof(:,p),...
                x_ups,x_fit_prof(:,p),...
                fit2DGauss.baselineEst,fit2DGauss.baselineEst,...
                [],imgData.blank,[],[]);      
            
        end
        
        eventParams.amplitude(2*i-1,1) = eP_2Dfit(1).amplitude;
        eventParams.amplitude(2*i,1) = eP_2Dfit(2).amplitude;
        
        eventParams.TTP(2*i-1,1) = eP_2Dfit(1).TTP;
        eventParams.TTP(2*i,1) = eP_2Dfit(2).TTP;
        
        eventParams.FDHM(2*i-1,1) = eP_2Dfit(1).FDHM;
        eventParams.FDHM(2*i,1) = eP_2Dfit(2).FDHM;
        
        eventParams.FWHM(2*i-1,1) = eP_2Dfit(1).FWHM;
        eventParams.FWHM(2*i,1) = eP_2Dfit(2).FWHM;
        
        eventParams.sparkMass(2*i-1,1) = eP_2Dfit(1).sparkMass;
        eventParams.sparkMass(2*i,1) = eP_2Dfit(2).sparkMass;

        pairedEventsFit(i) = {fit2DGauss};
        
        try
            % only if fit by EMG
            eventParams.tauD(2*i-1,1) = fit2DGauss.coefficientsOfFittedEvents{2,5};
            eventParams.tauD(2*i,1) = fit2DGauss.coefficientsOfFittedEvents{3,5};
            % normalize as (F-F0)/(F0-blank)
            D_fit_ups_norm1 = ...
                D_fit_ups_events{2,1}./(fit2DGauss.baselineEst - imgData.blank);
            D_fit_ups_norm2 = ...
                D_fit_ups_events{2,2}./(fit2DGauss.baselineEst - imgData.blank);
            % area under curve
            eventParams.AUC_2DFit(2*i-1,1) = ...
                sum(D_fit_ups_norm1(:).* mean(diff(t_ups)).*mean(diff(x_ups)));
            eventParams.AUC_2DFit(2*i,1) = ...
                sum(D_fit_ups_norm2(:).* mean(diff(t_ups)).*mean(diff(x_ups)));
            %eP_2Dfit.AUC_2DFit = trapz(x_ups,trapz(t_ups, D_fit_ups_norm, 2));  % area under curve
        catch
            eventParams.tauD(2*i-1,1) = nan;
            eventParams.tauD(2*i,1) = nan;
            eventParams.AUC_2DFit(2*i-1,1) = nan;
            eventParams.AUC_2DFit(2*i,1) = nan;
        end
        
    catch
        clearvars fit2DGauss
        eventParams.amplitude(2*i-1,1) = nan;
        eventParams.TTP(2*i-1,1) = nan;
        eventParams.FDHM(2*i-1,1) = nan;
        eventParams.FWHM(2*i-1,1) = nan;
        eventParams.sparkMass(2*i-1,1) = nan;
        eventParams.tauD(2*i-1,1) = nan;
        eventParams.AUC_2DFit(2*i-1,1) = nan;
        
        eventParams.amplitude(2*i,1) = nan;
        eventParams.TTP(2*i,1) = nan;
        eventParams.FDHM(2*i,1) = nan;
        eventParams.FWHM(2*i,1) = nan;
        eventParams.sparkMass(2*i,1) = nan;
        eventParams.tauD(2*i,1) = nan;
        eventParams.AUC_2DFit(2*i,1) = nan;
        
        pairedEventsFit(i) = {nan};
    end
    
    
    
    % show  individual calcium events
    if hObjs.check_showEventsFigs.Value
        
        %% create figure
        scF = 0.8;
        hf = figure('Tag','CaEventParam',...
            'Name',sprintf('posOfRoi %.2f um; pair #%d',posOfROI*pxSzX,i),...
            'Visible','on',...
            'OuterPosition',[(scrRes(3)-scrRes(3)*scF)/2 ...
            scrRes(2) ...
            scrRes(3)*scF ...
            scrRes(4)]);
        set(hf,'Units','Normalized')
        
        dx = 0.05;
        dy = 0.1;
        w_a = (1-5*dx)/4;
        h_a = 0.33;
        
        fontSzT = 14;
        fontSzL = 10;
        fontSzNum = 10;
        fontSzLegend = 12;
       
        % data to plot
        try
        firstEventFit = fit2DGauss.ups.individualEventsFits{2,1} + ...
            fit2DGauss.baselineEst;
        secondEventFit = fit2DGauss.ups.individualEventsFits{2,2} + ...
            fit2DGauss.baselineEst;
        catch
           keyboard 
        end
        %%%%% 2D gauss fit of pair of events, figure results %%%%%
        uicontrol('Style','text','Parent',hf,'Units','normalized',...
            'Position', [dx 1-0.5*dy 1-2*dx 0.5*dy],'FontUnits','normalized','FontSize',0.4,...
            'FontWeight','bold','HorizontalAlignment','left',...
            'String',sprintf('parameters from 2D exponentially modified gaussian fit (%s)',gauss2Dmodel));
        
        % show only relevant part of fit, where t profile is higher than 0.0001
        avrgTprof = mean( fit2DGauss.allEventsFit,1 ); 
        mShow = false(size(avrgTprof));
        try
            sShowMask = find( avrgTprof(1:find(avrgTprof > max(avrgTprof)/100,1,'first')) > 10E-5, ...
                1,'first');
            if isempty(sShowMask), sShowMask = 1; end
        catch
            sShowMask = 1;
        end
        
        % start of mask at least 20 ms before peak
        avrgTprofOfFirstEvent = mean(fit2DGauss.individualEventsFits{2,1},1);
        posOfMaxOfFirstEvent = find(avrgTprofOfFirstEvent == max(avrgTprofOfFirstEvent)); 
        if (posOfMaxOfFirstEvent - sShowMask)*pxSzT < 20 % ms
            sShowMask = round(posOfMaxOfFirstEvent - 20/pxSzT);
            if sShowMask < 1, sShowMask = 1; end  
        end
        eShowMask = find(avrgTprof > max(avrgTprof)/20,1,'last'); 
        mShow(sShowMask:eShowMask) = true;
        % upscaled
        mShowUps = false(size(t_ups));
        [~,sShowMaskUps] = min( abs(t_ups-t(sShowMask)) );
        [~,eShowMaskUps] = min( abs(t_ups-t(eShowMask)) );
        mShowUps(sShowMaskUps:eShowMaskUps) = true;
    
        ha2DGaussFit = axes(hf,'Units','normalized',...
            'Position',[dx 1-1.25*h_a-0.75*dy w_a*2.5 1.25*h_a]);
        % data
        mesh(T(:,mShow)-min(T(:,mShow),[],'all'), ...
            X(:,mShow), ...
            imgPair(:,mShow), ...
            'LineStyle','-',...
            'LineWidth',1,'FaceColor','none',...
            'EdgeColor','k','EdgeAlpha',0.4,...
            'FaceAlpha',0.2,'Parent',ha2DGaussFit)
        hold on
%         mesh(T_ups,X_ups,firstEventFit,...
%             'LineStyle','-','LineWidth',0.5,'FaceColor','b',...
%             'EdgeColor','none','EdgeAlpha',0.1,'FaceAlpha',0.4)
%         mesh(T_ups,X_ups,secondEventFit,'LineStyle','-',...
%             'LineWidth',0.5,'FaceColor','none',...
%             'EdgeColor','g','EdgeAlpha',0.1,'FaceAlpha',0.4)

        % fit
        mesh(T_ups(:,mShowUps)-min(T_ups(:,mShowUps),[],'all'), ...
            X_ups(:,mShowUps),...
            fit2DGauss.ups.wholeFit(:,mShowUps),...
            'LineStyle','-',...
            'LineWidth',0.5,'FaceColor','interp',...
            'EdgeColor','none','EdgeAlpha',0.1,'FaceAlpha',0.4)
        colormap(ha2DGaussFit,jet(256))

        line(t_ups(mShowUps)-min(t_ups(mShowUps)),...
            X_ups(rFitMaxUps(1),mShowUps),t_fit_prof(mShowUps,1),...
            'Parent',ha2DGaussFit,'LineWidth',3,'Color','b')
        line(T_ups(:,cFitMaxUps(1))-min(t_ups(mShowUps)),x_ups,x_fit_prof(:,1),...
            'Parent',ha2DGaussFit,'LineWidth',3,'Color','b')
        line(t_ups(mShowUps)-min(t_ups(mShowUps)),X_ups(rFitMaxUps(2),mShowUps),t_fit_prof(mShowUps,2),...
            'Parent',ha2DGaussFit,'LineWidth',3,'Color','k')
        line(T_ups(:,cFitMaxUps(2))-min(t_ups(mShowUps)),x_ups,x_fit_prof(:,2),...
            'Parent',ha2DGaussFit,'LineWidth',3,'Color','k')
        
        % add TPP as patches
        % make sure TPP are 1 ms long make eps or pdf to be able to manipulate in 
        % position of TPP pulses in x axis
        
        % put position of 2PP in the position of max of first event
        TPP_x1 = posOfTTP_inPairEventsArea - imgData.d_TP_laser/2;
        TPP_x2 = posOfTTP_inPairEventsArea + imgData.d_TP_laser/2;
        
        TPP_x1 = mean(X_ups(rFitMaxUps(1),mShowUps)) - imgData.d_TP_laser/2;
        TPP_x2 = mean(X_ups(rFitMaxUps(1),mShowUps)) + imgData.d_TP_laser/2;
        
        sPP_pair = sPP_pair - min(t(mShow));
        ePP_pair = ePP_pair - min(t(mShow));
        ePP_pair = sPP_pair + 1;
        vertices_2PP = [ ...
            sPP_pair(1) TPP_x1 ha2DGaussFit.ZLim(1); ...
            ePP_pair(1) TPP_x1 ha2DGaussFit.ZLim(1); ...
            ePP_pair(1) TPP_x2 ha2DGaussFit.ZLim(1); ...
            sPP_pair(1) TPP_x2 ha2DGaussFit.ZLim(1); ...
            sPP_pair(1) TPP_x2 ha2DGaussFit.ZLim(2); ...
            ePP_pair(1) TPP_x2 ha2DGaussFit.ZLim(2); ...
            ePP_pair(1) TPP_x1 ha2DGaussFit.ZLim(2); ...
            sPP_pair(1) TPP_x1 ha2DGaussFit.ZLim(2); ...
            sPP_pair(2) TPP_x1 ha2DGaussFit.ZLim(1); ...
            ePP_pair(2) TPP_x1 ha2DGaussFit.ZLim(1); ...
            ePP_pair(2) TPP_x2 ha2DGaussFit.ZLim(1); ...
            sPP_pair(2) TPP_x2 ha2DGaussFit.ZLim(1); ...
            sPP_pair(2) TPP_x2 ha2DGaussFit.ZLim(2); ...
            ePP_pair(2) TPP_x2 ha2DGaussFit.ZLim(2); ...
            ePP_pair(2) TPP_x1 ha2DGaussFit.ZLim(2); ...
            sPP_pair(2) TPP_x1 ha2DGaussFit.ZLim(2); ...
            ];
        %3D
        faces_2PP_1 = [ ...
            1 2 3 4; ...
            2 3 6 7; ...
            5 6 7 8; ...
            1 2 7 8; ...
            1 4 5 8; ...
            3 4 5 6 ];
        faces_2PP_2 = faces_2PP_1 + 8;
        faces_2PP = [faces_2PP_1;faces_2PP_2];
%         % 2D
%         faces_2PP = [1 2 3 4;9 10 11 12];
        
        patch(ha2DGaussFit,'Faces',faces_2PP,'Vertices',vertices_2PP, ...
            'FaceColor','r','FaceAlpha',0.25,'EdgeColor','none')
     
        view(ha2DGaussFit,-5,15)
        set(ha2DGaussFit,'YDir','normal','FontSize',fontSzNum)
        title(ha2DGaussFit,'ff')
        title(ha2DGaussFit,...
            sprintf('posOfRoi %.2f \\mum; pair of calcium events (#%d); 2D gauss fit ',posOfROI*pxSzX,i),...
            'FontSize',fontSzT)
        xlabel(ha2DGaussFit,'t (ms)','FontSize',fontSzL)
        ylabel(ha2DGaussFit,'x (\mum)','FontSize',fontSzL)
        zlabel(hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
        set(ha2DGaussFit, ...
            'XLim',[min(t(mShow)-min(t(mShow))), max(t(mShow)-min(t(mShow)))], ...
            'YLim',[min(x) max(x)], ...
            'ZLim',[min([imgPair(:);fit2DGauss.ups.wholeFit(:)])*0.99 ... 
            max([imgPair(:);fit2DGauss.ups.wholeFit(:)])*1.01])
        % remove grids
        set(ha2DGaussFit, 'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off')
        
        set(ha2DGaussFit, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on')

          keyboard
        % visibility
      
        ha2DGaussFit.Children(1).Visible = 'off'
        ha2DGaussFit.Children(6).Visible = 'on'
        ha2DGaussFit.Children(7).Visible = 'off'
        ha2DGaussFit.Children(2).Visible = 'on'
        ha2DGaussFit.Children(3).Visible = 'on'
        ha2DGaussFit.Children(4).Visible = 'on'
        ha2DGaussFit.Children(5).Visible = 'on'
        
        % axes for x,t profiles for individual 2d gauss fits in pair
        ha_e1_xP = axes(hf,'Units','normalized',...
            'Position',[dx 1-2*h_a-2.5*dy w_a h_a]);
        line(x_ups,x_fit_prof(:,1),'Parent',ha_e1_xP,'LineWidth',2,...
            'Color','b','Display','fit')
        set(ha_e1_xP,'FontSize',fontSzNum)
        title(ha_e1_xP,sprintf('1. calcium event in pair #%d; x profile',i),...
            'FontSize',fontSzT)
        ylabel(ha_e1_xP,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
        xlabel(ha_e1_xP,'x (\mum)','FontSize',fontSzL)
        xlim(ha_e1_xP,[min(x_ups) max(x_ups)])
        ylim(ha_e1_xP,[min(x_fit_prof(:,1))*0.95 max(x_fit_prof(:,1))*1.05])
        line([eP_2Dfit(1).half_max_x_1 eP_2Dfit(1).half_max_x_2],...
            [eP_2Dfit(1).half_max_x eP_2Dfit(1).half_max_x],'Parent',ha_e1_xP,...
            'LineWidth',2,'Color','r','Display','FWHM')
        hl_2 = legend(ha_e1_xP,'show');
        hl_2.Location = 'northwest';
        hl_2.FontSize = fontSzLegend;
  
        % 
        ha_e1_tP = axes(hf,'Units','normalized',...
            'Position',[2*dx+w_a 1-2*h_a-2.5*dy w_a h_a]);
        line(t_ups(mShowUps),t_fit_prof(mShowUps,1),'Parent',ha_e1_tP,'LineWidth',2,...
            'Color','b','Display','fit')
        set(ha_e1_tP,'FontSize',fontSzNum)
        title(ha_e1_tP,sprintf('1. calcium event in pair #%d; t profile',i),...
            'FontSize',fontSzT)
        xlabel(ha_e1_tP,'t (ms)','FontSize',fontSzL)
        ylabel(ha_e1_tP,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
        xlim(ha_e1_tP,[min(t_ups(mShowUps)) max(t_ups(mShowUps))])
        ylim(ha_e1_tP,[min(t_fit_prof(:,1))*0.95 max(t_fit_prof(:,1))*1.05])
        line([eP_2Dfit(1).half_max_t_1 eP_2Dfit(1).half_max_t_2],...
            [eP_2Dfit(1).half_max_t eP_2Dfit(1).half_max_t],'Parent',ha_e1_tP,...
            'LineWidth',2,'Color','r','Display','FDHM')
        line([eP_2Dfit(1).t_max eP_2Dfit(1).t_max],...
            [eP_2Dfit(1).bs_t eP_2Dfit(1).v_max],'Parent',ha_e1_tP,...
            'LineWidth',2,'Color','g','Display','amplitude')
        line([eP_2Dfit(1).t0 eP_2Dfit(1).t_max],...
            [eP_2Dfit(1).bs_t eP_2Dfit(1).bs_t],'Parent',ha_e1_tP,...
            'LineWidth',2,'Color','m','Display','TTP')
        hl_3 = legend(ha_e1_tP,'show');
        hl_3.Location = 'best';
        hl_3.FontSize = fontSzLegend;
        
        %
        ha_e2_xP = axes(hf,'Units','normalized',...
            'Position',[3*dx+2*w_a 1-2*h_a-2.5*dy w_a h_a]);
        line(x_ups,x_fit_prof(:,2),'Parent',ha_e2_xP,'LineWidth',2,...
            'Color','k','Display','fit')
        set(ha_e2_xP,'FontSize',fontSzNum)
        title(ha_e2_xP,sprintf('2. calcium event in pair #%d; x profile',i),...
            'FontSize',fontSzT)
        ylabel(ha_e2_xP,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
        xlabel(ha_e2_xP,'x (\mum)','FontSize',fontSzL)
        xlim(ha_e2_xP,[min(x_ups) max(x_ups)])
        try
            ylim(ha_e2_xP,[min(x_fit_prof(:,2))*0.95 max(x_fit_prof(:,2))*1.05])
        catch 
            keyboard
            ylim(ha_e2_xP,[min(x_fit_prof(:,2))*0.95-1 max(x_fit_prof(:,2))*1.05+1])
        end
        line([eP_2Dfit(2).half_max_x_1 eP_2Dfit(2).half_max_x_2],...
            [eP_2Dfit(2).half_max_x eP_2Dfit(2).half_max_x],'Parent',ha_e2_xP,...
            'LineWidth',2,'Color','r','Display','FWHM')
        hl_4 = legend(ha_e2_xP,'show');
        hl_4.Location = 'northwest';
        hl_4.FontSize = fontSzLegend;
  
        %
        ha_e2_tP = axes(hf,'Units','normalized',...
            'Position',[4*dx+3*w_a 1-2*h_a-2.5*dy w_a h_a]);
        line(t_ups(mShowUps),t_fit_prof(mShowUps,2),'Parent',ha_e2_tP,'LineWidth',2,...
            'Color','k','Display','fit')
        set(ha_e2_tP,'FontSize',fontSzNum)
        title(ha_e2_tP,sprintf('2. calcium event in pair #%d; t profile',i),...
            'FontSize',fontSzT)
        xlabel(ha_e2_tP,'t (ms)','FontSize',fontSzL)
        ylabel(ha_e2_tP,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
        xlim(ha_e2_tP,[min(t_ups(mShowUps)) max(t_ups(mShowUps))])
        ylim(ha_e2_tP,[min(t_fit_prof(:,2))*0.95 max(t_fit_prof(:,2))*1.05])
        line([eP_2Dfit(2).half_max_t_1 eP_2Dfit(2).half_max_t_2],...
            [eP_2Dfit(2).half_max_t eP_2Dfit(2).half_max_t],'Parent',ha_e2_tP,...
            'LineWidth',2,'Color','r','Display','FDHM')
        line([eP_2Dfit(2).t_max eP_2Dfit(2).t_max],...
            [eP_2Dfit(2).bs_t eP_2Dfit(2).v_max],'Parent',ha_e2_tP,...
            'LineWidth',2,'Color','g','Display','amplitude')
        line([eP_2Dfit(2).t0 eP_2Dfit(2).t_max],...
            [eP_2Dfit(2).bs_t eP_2Dfit(2).bs_t],'Parent',ha_e2_tP,...
            'LineWidth',2,'Color','m','Display','TTP')
        hl_5 = legend(ha_e2_tP,'show');
        hl_5.Location = 'best';
        hl_5.FontSize = fontSzLegend;
        
        % text
        txtStr = string(sprintf('parameters of calcium events from pair #%d:',i));

        paramsNames = [
            string(['amplitude',' (',char(916),'F/F_0)']);
            "TTP (ms)";
            "FDHM (ms)";
            string(['FWHM',' (',char(181),'m)']);
            string(['sparkMass',' (',char(916),'F/F_0*',char(181),'m^3) ']);
            string(['volume under fit',' (',char(916),'F/F_0*',char(181),'m*ms) '])
            ];
        paramsNames = pad(paramsNames,' ');
        
        paramsNames = [
            pad("parameters:",strlength(paramsNames(1)));
            pad("-",strlength(paramsNames(1)),'-'); 
            paramsNames
            ];
        % add spaces to compensate for super and subscripts
        paramsNames(3) = join([paramsNames(3),""]," ");
        paramsNames(7) = join([paramsNames(7),""],"  ");
        paramsNames(8) = join([paramsNames(8),""]," ");
        
        paramsValEvent1 = pad([
            "   1. event   ";
            "--------------";
            string(sprintf('%0.2f',eP_2Dfit(1).amplitude));
            string(sprintf('%0.2f',eP_2Dfit(1).TTP));
            string(sprintf('%0.2f',eP_2Dfit(1).FDHM));
            string(sprintf('%0.2f',eP_2Dfit(1).FWHM));
            string(sprintf('%0.2f',eP_2Dfit(1).sparkMass));
            string(sprintf('%0.2f',eventParams.AUC_2DFit(2*i-1,1)))
            ],'both');

        paramsValEvent2 = pad([
            "   2. event   ";
            "--------------";
            string(sprintf('%0.2f',eP_2Dfit(2).amplitude));
            string(sprintf('%0.2f',eP_2Dfit(2).TTP));
            string(sprintf('%0.2f',eP_2Dfit(2).FDHM));
            string(sprintf('%0.2f',eP_2Dfit(2).FWHM));
            string(sprintf('%0.2f',eP_2Dfit(2).sparkMass));
            string(sprintf('%0.2f',eventParams.AUC_2DFit(2*i,1)))
            ],'both');
       
        params = join( [ strings(size(paramsNames)),paramsNames,...
            paramsValEvent1,paramsValEvent2,strings(size(paramsNames)) ],'|' );
        
        txtStr = [ 
            pad(txtStr,strlength(params(1))); 
            pad("-",strlength(params(1)),'-');
            params ];
        
        % invisible axes for text
        ha_txt = axes(hf,'Units','normalized',...
            'Position',[ha2DGaussFit.Position(1)+ha2DGaussFit.Position(3)+dx ...
            ha2DGaussFit.Position(2)...
            1-ha2DGaussFit.Position(1)-ha2DGaussFit.Position(3)-dx ...
            ha2DGaussFit.Position(4)]);
        set(ha_txt,'Visible','off')
        text(min(get(ha_txt,'XLim')),max(get(ha_txt,'YLim')),txtStr,...
        'Parent',ha_txt,'FontUnits','points','FontSize',12,...
        'VerticalAlignment','top','FontName','FixedWidth')
        

    end
    
    keyboard % try to save it so it can be broken into parts in adobe illustator or canvas
    
    clearvars D t x t_prof x_prof t_ups x_ups t_prof_ups x_prof_ups ... 
        y_t1 y_t2 D1 y_x1 y_x2
    
    clearvars rows_e cols_e eP_profs eP_2Dfit out_t_prof out_x_prof p0 ...
        t_spark_prof x_spark_prof fit2DGauss
    
    % Report current estimate in the waitbar's message field
    waitbar(i/N_pairs,hw,sprintf('%d %%',round(100*i/N_pairs)))
    
end

% delete waitbar
delete(hw), clearvars hw

end


function allFits = findDetectedEventsParamsLocalFun(...
    img,statEvents,startOfEvent,endOfEvent,prevFitCoef,...
    mainFig,calcMethod,gauss2Dmodel)

% img = not normalized data, raw data
% statEvents = statistic from event detection
% startOfSpark,endOfSpark,prevFitCoef = from previous fitting in spark recovery analysis 

% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;

% filter raw image
img = imgFiltering(img,pxSzT,pxSzX);

h_d = str2double(get(hObjs.h_edit_averageWidth,'String'));

n_px_t = ceil(h_d/pxSzX);
if mod(n_px_t,2)==0
    n_px_t = n_px_t - 1;
end

n_px_x = ceil(5/pxSzT);
if mod(n_px_x,2)==0
    n_px_x = n_px_x - 1;
end

% out = eventBoundaries(img,statEvents,mainFig)
% % progress bar 
% progressbar('calculating detected events parameters')

for i = 1:length(statEvents)
    
    [imgE,imgEs,imgE_m,maxOfEventPos,maxCrossProfs,t0,bs] = ...
        eventImgAndMaxCrossingProfiles( ...
        statEvents(i),startOfEvent(i),endOfEvent(i),...
        pxSzT,pxSzX,n_px_t,n_px_x,img,prevFitCoef(i,:),mainFig);

    % choose method how to calculate spark parameters
    switch calcMethod
        
        case '2DGauss'
            % calculate parameters, fit spark data with 2D exponentially modified gaussian
            % create time and spatial grids
            [T,X] = meshgrid(maxCrossProfs.t,maxCrossProfs.x);
            [T_ups,X_ups] = meshgrid(maxCrossProfs.t_ups,maxCrossProfs.x_ups);

            try
                fit2DGauss = fit2DGaussian( ...
                    T,X,T_ups,X_ups,imgE,imgE_m,maxOfEventPos,gauss2Dmodel);
                fit2DGauss.dataEventMask = imgE_m;
            catch
                
                fit2DGauss = [];
            end
            allFits(i,1) = {fit2DGauss};
            
        case 'peakXTProfiles'
            
            try
                % get parameters of spark from profiles without 2D fitting
                % try fitting
                % starting points for chosen fit function
                % coefEst = paramEstimEMG(t_prof,t,t0)
                % 'EMG' : coeff = [A,m,sd,tau,F0]
                p0 = [maxCrossProfs.peakData_tProf.val ...
                    maxCrossProfs.peakData_tProf.pos*pxSzT 5 15 bs];
                out_t_prof = fitOneWholeEvent( ...
                    p0,maxCrossProfs.t,maxCrossProfs.t_ups,...
                    maxCrossProfs.t_event_prof,'no','EMG',...
                    maxCrossProfs.t_event_prof_m);
                
                % fit exp modified gaussian or spline
                % [A,m,sd,tau,F0]
                xProf_m = imgE_m(:,maxOfEventPos(2));
                x_event_prof_m = maxCrossProfs.x_event_prof;
                x_event_prof_m(~xProf_m) = nan;
                [v_max_x,p_max_x] = max(x_event_prof_m);
               
                out_x_prof = fitOneWholeEvent( ...
                    [v_max_x p_max_x*pxSzX-2 1 3 bs], ...
                    maxCrossProfs.x,maxCrossProfs.x_ups,...
                    maxCrossProfs.x_event_prof,'no','EMG',xProf_m);
 
%                 figure
%                 plot(out_t_prof.x,out_t_prof.y,'ok')
%                 hold on
%                 plot(out_t_prof.x,out_t_prof.yFit,'r')
                     
            catch
                % if cannot get profiles for some reason
                out_t_prof = [];
                out_x_prof = [];
                
            end
            allFits(i,:) = {out_t_prof,out_x_prof};
    end

    clearvars imgE t x t_prof x_prof t_ups x_ups t_prof_ups x_prof_ups y_t1 y_t2 D1 y_x1 y_x2
    
    clearvars rows_e cols_e eP_profs eP_2Dfit out_t_prof out_x_prof p0 t_spark_prof x_spark_prof

end

end


