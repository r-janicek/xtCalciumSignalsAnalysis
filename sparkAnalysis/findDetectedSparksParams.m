function eventParams = findDetectedSparksParams(...
    img, statSparks, mainFig ,calcMethod, sparksRec, posOfROI, ...
    startOfSpark, endOfSpark, prevFitCoeffs, normImgFlag)

% img = filtered data or filtered normalized data
% statSparks = statistic from sparks detection
% startOfSpark, endOfSpark, prevFitCoeffs = previous fit with fitSparkRise
% function

% posOfROI = from spark recovery analysis
if isempty(posOfROI)
    posOfROI = nan;
end

% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;
scrRes = get(0,'ScreenSize');

% t-profile = average from 1 um width area
% x-profile = average from 5 ms width area
n_px_t = ceil(1/pxSzX);
if mod(n_px_t,2)==0
    n_px_t = n_px_t - 1;
end

n_px_x = ceil(5/pxSzT);
if mod(n_px_x,2)==0
    n_px_x = n_px_x - 1;
end

% out = eventBoundaries(img,statEvents,mainFig)
% change pointer to watch
set(mainFig,'Pointer','watch')
drawnow

% check if there was some sparks detected
if ~isempty(statSparks)
    % add waitbar
    hw = waitbar(0,'0 %','Name','analyzing detected sparks...',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(hw,'canceling',0)
    
    if isfield(getappdata(mainFig),'batchProcessing')
        if getappdata(mainFig,'batchProcessing')
            hw.Visible = 'off';
        end
    end
    
    for i = 1:length(statSparks)
        % get number of spark ROI rectangle
        switch getappdata(mainFig, 'analysisType')
            case 'spark detection'
                sparkROInum = str2double(sparksRec(i).Tag);
            case 'spark recovery ryanodine'
                sparkROInum = i;
        end
        if isempty(startOfSpark)
            sSpark = [];
            eSpark = [];
            prevSparkFitCoeffs = [];
        else
            sSpark = startOfSpark(i);
            eSpark = endOfSpark(i);
            prevSparkFitCoeffs = prevFitCoeffs(i,:);
        end
        % get image of event and max crossing profiles
        [imgE, imgEs, imgE_m, maxOfEventPos, ... 
            maxCrossProfs, t0, bs, eventROIstart_t] = ...
            eventImgAndMaxCrossingProfiles( ...
            statSparks(i), sSpark, eSpark,...
            pxSzT, pxSzX, n_px_t, n_px_x, img, ...
            prevSparkFitCoeffs, mainFig);
        
        t = maxCrossProfs.t;
        t_ups = maxCrossProfs.t_ups;
        t_event_prof = maxCrossProfs.t_event_prof;
        t_event_profS = maxCrossProfs.t_event_prof_S;
        
        peakData_tProf = maxCrossProfs.peakData_tProf;
        
        x = maxCrossProfs.x;
        x_ups = maxCrossProfs.x_ups;
        x_event_prof = maxCrossProfs.x_event_prof;
        x_event_profS = maxCrossProfs.x_event_prof_S;
        
        r_m = maxOfEventPos(1);
        c_m = maxOfEventPos(2);
        % masks of profiles
        tProf_m = maxCrossProfs.t_event_prof_m;
        xProf_m = imgE_m(:,c_m);
            
        % choose method how to calculate spark parameters
        calcMethodPrev = calcMethod;
        switch calcMethod
            
            case '2DGauss'
                % calculate parameters, fit spark data with 2D exponentially modified gaussian
                try
                    try
                        % create time and spatial grids
                        [T,X] = meshgrid(t,x);
                        [T_ups,X_ups] = meshgrid(t_ups,x_ups);
                        
                        gaussModel = '2D_EMGxt';
                        fit2DGauss = fit2DGaussian(...
                            T,X,T_ups,X_ups,imgE,imgE_m,maxOfEventPos,gaussModel);
                        
                        imgE_fit = fit2DGauss.dataEventFit;
                        imgE_fit_ups = fit2DGauss.dataEventFit_ups;
                        
                        % get t and x profiles crossing max of fit
                        [rFitMax,cFitMax] = find(imgE_fit == max(imgE_fit(:)));
                        t_event_prof = imgE(rFitMax,:); t_event_prof = t_event_prof(:);
                        x_event_prof = imgE(:,cFitMax); x_event_prof = x_event_prof(:);
                        [rFitMaxUps,cFitMaxUps] = find(imgE_fit_ups == max(imgE_fit_ups(:)));
                        t_prof_fit = imgE_fit_ups(rFitMaxUps,:);
                        x_prof_fit = imgE_fit_ups(:,cFitMaxUps);
                        t_prof_fit = t_prof_fit(:);
                        x_prof_fit = x_prof_fit(:);
                        
                        [~,peakData_ups.pos] = ...
                            min( abs( t_ups - t(maxCrossProfs.peakData_tProf.pos) ) );
                        peakData_ups.val = t_prof_fit(peakData_ups.pos);
                        
                        eP_2Dfit = getParametersOfEventProfile( ...
                            t_ups, t_prof_fit,...
                            x_ups, x_prof_fit,...
                            t_prof_fit(1), t_prof_fit(1), [],...
                            imgData.blank, peakData_ups, [], normImgFlag);
                         
                    catch
                        % get params from image of event
                        calcMethodPrev = '2DGauss'; 
                        calcMethod = 'estimate from event img';
                        
                        eP_2Dfit = getParametersOfEventFromRawImage(...
                            imgE, imgEs, imgE_m, ...
                            [], [], r_m, c_m, ...
                            t_event_profS, x_event_profS, t, x, ...
                            imgData, normImgFlag);

                    end
                    
                    eventParams.amplitude(i,1) = eP_2Dfit.amplitude;
                    eventParams.t0(i,1) = eP_2Dfit.t0 + eventROIstart_t;
                    eventParams.TTP(i,1) = eP_2Dfit.TTP;
                    eventParams.FDHM(i,1) = eP_2Dfit.FDHM;
                    eventParams.FWHM(i,1) = eP_2Dfit.FWHM;
                    eventParams.sparkMass(i,1) = eP_2Dfit.sparkMass;
                    
                    try
                        % only if fit by EMG
                        eventParams.tauD(i,1) = fit2DGauss.coef2DGauss(4);
                        % (F-F0)/(F0-blank)
                        imgE_fit_ups_norm = ...
                            (imgE_fit_ups - fit2DGauss.coef2DGauss(end))/ ...
                            (fit2DGauss.coef2DGauss(end) - imgData.blank);
                        % area under curve
                        eventParams.AUC_2DFit(i,1) = ...
                            sum(imgE_fit_ups_norm(:).* mean(diff(t_ups)).*mean(diff(x_ups)));
                        %eP_2Dfit.AUC_2DFit = trapz(x_ups,trapz(t_ups, D_fit_ups_norm, 2));
                    catch
                        eventParams.AUC_2DFit(i,1) = nan;
                    end
                     
                catch
                    eventParams.amplitude(i,1) = nan;
                    eventParams.t0(i,1) = nan;
                    eventParams.TTP(i,1) = nan;
                    eventParams.FDHM(i,1) = nan;
                    eventParams.FWHM(i,1) = nan;
                    eventParams.sparkMass(i,1) = nan;
                    eventParams.tauD(i,1) = nan;
                    eventParams.AUC_2DFit(i,1) = nan;
                end
                  
            case 'peakXTProfiles'
                      
                try
                    % get parameters of spark from profiles without 2D fitting            
                    % do fits of profiles
                    try
                        % try fitting
                        % '1expR1expD' piecewise: coefficients: [t0, F01, tauR, A, t1, tauD, F02]
                        % 'spline' piecewise: coefficients: [t0, F01]
                        % modified gaussian = 'EMG', coefficients:[A,m,sd,tau,F0]
                        % p0 = [t0 bs tauR v_max_t p_max_t*pxSzT 15 mean(t_spark_prof(end-5:end))];
                        % p0 = [v_max_t p_max_t*pxSzT p_max_t*pxSzT-t0 15 bs];
                        % initial parameters for time profile fit
                        p0 = [t0 bs];
                        out_t_prof = fitOneWholeEvent( ...
                            p0, t, t_ups, t_event_prof, ...
                            'yes', 'spline', tProf_m);
                        
                        bs_t_prof = out_t_prof.bs;
                        t0_t_prof = out_t_prof.t0;
                        tauD = out_t_prof.tauD;
                        
                        % fit gaussian or exp modified gaussian or spline
                        % to spatial profile
                        % 'Gauss' : [F0,A,w,xc]
                        % 'EMG' : [A,m,sd,tau,F0]
                        x_event_prof_m = x_event_prof;
                        x_event_prof_m(~xProf_m) = nan;
                        [v_max_x, p_max_x] = max(x_event_prof_m);
                        
                        % estimation of the width in half max
                        w = sum((x_event_prof_m-bs-(v_max_x-bs)/2)>=0)*pxSzX;
                        % fit
                        out_x_prof = fitOneWholeEvent(...
                            [bs v_max_x-bs w p_max_x*pxSzX-pxSzX],...
                            x, x_ups, x_event_prof, 'no', 'Gauss', xProf_m);
                        
                        % baseline of x profile
                        bs_x_prof = out_x_prof.bs;
                        
                        % upscale mask
                        % start and end of event's mask
                        tProf_m_ups = false(size(t_ups));
                        sE = find(tProf_m,1,'first');
                        eE = find(tProf_m,1,'last');
                        [~, sE_ups] = min( abs( t_ups - t(sE) ) );
                        [~, eE_ups] = min( abs( t_ups - t(eE) ) );
                        tProf_m_ups(sE_ups:eE_ups) = true;
                        % peak position and value
                        [~, loc_ups] = findpeaks(...
                            out_t_prof.yFit_ups(tProf_m_ups), ...
                            t_ups(tProf_m_ups), 'NPeaks',1, ...
                            'SortStr','descend');
                        [~, peakData_ups.pos] = ...
                            min( abs( t_ups - loc_ups ) );
                        peakData_ups.val = ...
                            out_t_prof.yFit_ups(peakData_ups.pos);

                        eP_profs = getParametersOfEventProfile( ...
                            t_ups, out_t_prof.yFit_ups,...
                            x_ups, out_x_prof.yFit_ups,...
                            bs_t_prof, bs_x_prof, t0_t_prof,...
                            imgData.blank, peakData_ups, tProf_m_ups, ...
                            normImgFlag);
                       
                    catch
                        % get params from image of event
                        calcMethodPrev = 'peakXTProfiles'; 
                        calcMethod = 'estimate from event img';
                        
                        eP_profs = getParametersOfEventFromRawImage(...
                            imgE, imgEs, imgE_m, ...
                            [], [], r_m, c_m, ...
                            t_event_profS, x_event_profS, t, x, ...
                            imgData, normImgFlag);
                    end
                    
                    eventParams.amplitude(i,1) = eP_profs.amplitude;
                    eventParams.t0(i,1) = eP_profs.t0 + eventROIstart_t;
                    eventParams.TTP(i,1) = eP_profs.TTP;
                    eventParams.FDHM(i,1) = eP_profs.FDHM;
                    eventParams.FWHM(i,1) = eP_profs.FWHM;
                    eventParams.sparkMass(i,1) = eP_profs.sparkMass;
                    eventParams.tauD(i,1) = tauD;
                    eventParams.AUC_2DFit(i,1) = nan;
                      
                catch
                    % if cannot get profiles for some reason
                    eventParams.amplitude(i,1) = nan;
                    eventParams.t0(i,1) = nan;
                    eventParams.TTP(i,1) = nan;
                    eventParams.FDHM(i,1) = nan;
                    eventParams.FWHM(i,1) = nan;
                    eventParams.sparkMass(i,1) = nan;
                    eventParams.tauD(i,1) = nan;
                    eventParams.AUC_2DFit(i,1) = nan;
                    
                end
        
                
            case 'estimate from event img'
                try
                    % get parameters of spark from profiles without 2D fitting
                    % use coeficient from previous fitting (t0, tauR, A, baseline)
                    % get params from image of event
                    if ~isempty(t0)
                        eParamsEst = getParametersOfEventProfile( ...
                            t, t_event_profS,...
                            x, x_event_profS,...
                            bs, bs, t0,...
                            imgData.blank, ...
                            peakData_tProf, tProf_m, normImgFlag);
                        % estimate tauD
                        try
                            [~, max_pos_s] = ...
                                max(t_event_profS(maxCrossProfs.t_event_prof_m));
                            max_pos_s = max_pos_s - 1 + ... 
                                find(maxCrossProfs.t_event_prof_m, 1, 'first');
                            profD = t_event_profS(max_pos_s+1:end);
                            try
                                % last 3 points
                                bs_end = mean(profD(end-2:end)); 
                            catch
                                % last point
                                bs_end = profD(end); 
                            end
                            [~, p_tauD] = min( abs( ...
                                profD - (0.33*(t_event_profS(max_pos_s)-bs_end)+bs_end) ) );
                            tauD = p_tauD * imgData.pxSzT;
                        catch
                            tauD = nan;
                        end
                        % eParamsEst = ...
                        %     getParametersOfEventFromRawImage(...
                        %     imgE, imgEs, imgE_m, ...
                        %     bs, t0, r_m, c_m, ...
                        %     t_event_profS, x_event_profS, t, x, imgData);
                    else
                        eParamsEst = ...
                            getParametersOfEventFromRawImage(...
                            imgE, imgEs, imgE_m, ...
                            [], [], r_m, c_m, ...
                            t_event_profS, x_event_profS, t, x, ...
                            imgData, normImgFlag);
                    end
                    
                    eventParams.amplitude(i,1) = eParamsEst.amplitude;
                    eventParams.t0(i,1) = eParamsEst.t0 + eventROIstart_t;
                    eventParams.TTP(i,1) = eParamsEst.TTP;
                    eventParams.FDHM(i,1) = eParamsEst.FDHM;
                    eventParams.FWHM(i,1) = eParamsEst.FWHM;
                    eventParams.sparkMass(i,1) = eParamsEst.sparkMass;
                    eventParams.tauD(i,1) = tauD;
                    eventParams.AUC_2DFit(i,1) = nan;
                    
                catch
                    % if cannot get estimate of params for some reason
                    eventParams.amplitude(i,1) = nan;
                    eventParams.t0(i,1) = nan;
                    eventParams.TTP(i,1) = nan;
                    eventParams.FDHM(i,1) = nan;
                    eventParams.FWHM(i,1) = nan;
                    eventParams.sparkMass(i,1) = nan;
                    eventParams.tauD(i,1) = nan;
                    eventParams.AUC_2DFit(i,1) = nan;
                end
                
        end
        
        eventParams.sparkROINum(i,1) = sparkROInum;
        
        % show  individual calcium events
        if hObjs.check_showEventsFigs.Value
            % setup titles of figures and axes
            switch getappdata(mainFig, 'analysisType')
                case 'spark detection'
                    fig_title_txt = ...
                        sprintf('spark from ROI #: %d',sparkROInum);
                    ax_profX_title_txt = ...
                        sprintf('spark from ROI #: %d -- x profile',sparkROInum);
                    ax_profT_title_txt = ...
                        sprintf('spark from ROI #: %d -- t profile',sparkROInum);
                    txt_b = {sprintf('parameters of spark from ROI #%d:', ...
                        sparkROInum)};
                case 'spark recovery ryanodine'
                    fig_title_txt = ...
                        sprintf('posOfRoi %.2f \\mum; spark #%d', ...
                        posOfROI*pxSzX, sparkROInum);
                    ax_profX_title_txt = ...
                        sprintf('posOfRoi %.2f \\mum; spark #%d -- x profile', ...
                        posOfROI*pxSzX, sparkROInum);
                    ax_profT_title_txt = ...
                        sprintf('posOfRoi %.2f \\mum; spark #%d -- t profile', ...
                        posOfROI*pxSzX, sparkROInum);
                    txt_b = {[sprintf('posOfRoi %.2f ',posOfROI*pxSzX), ...
                        char(181),'m -- ', sprintf('parameters of spark #%d:', ...
                        sparkROInum)]};
            end
            % create figure
            scF = 0.8;
            hf = figure('Tag','CaEventParam',...
                'Name',fig_title_txt,...
                'Visible','on',...
                'OuterPosition',[scrRes(3)-scrRes(3)*scF ...
                scrRes(4)*0.4 ...
                scrRes(3)*scF ...
                scrRes(4)*0.6]);
            set(hf,'Units','Normalized')
            
            h_a = 0.55;
            w_a = 0.225;
            dx = 0.05;
            dy = 0.07;
            
            fontSzT = 18;
            fontSzL = 16;
            fontSzNum = 14;
            fontSzLegend = 14;

            switch calcMethod
                case '2DGauss'
                    txt_a = ['parameters from', ... 
                        '2D exponentially modified gaussian fit', ... 
                        sprintf('(%s)',gaussModel)];
                    x_prof_fit_txt = 'fit';
                    half_max_x_1 = eP_2Dfit.half_max_x_1;
                    half_max_x_2 = eP_2Dfit.half_max_x_2;
                    half_max_x = eP_2Dfit.half_max_x;
                    t_prof_fit_txt = 'fit';
                    half_max_t_1 = eP_2Dfit.half_max_t_1;
                    half_max_t_2 = eP_2Dfit.half_max_t_2;
                    half_max_t = eP_2Dfit.half_max_t;
                    t_max = eP_2Dfit.t_max;
                    bs_t = eP_2Dfit.bs_t;
                    v_max = eP_2Dfit.v_max;
                    t0 = eP_2Dfit.t0;
    
                case 'peakXTProfiles'
                    txt_a = ['parameters from profiles',...
                        '(1 um and 5 ms wide) crossing peak of event'];
                    x_prof_fit = out_x_prof.yFit_ups;
                    x_prof_fit_txt = ...
                        sprintf('fit with %s', out_x_prof.fitModel);
                    half_max_x_1 = eP_profs.half_max_x_1;
                    half_max_x_2 = eP_profs.half_max_x_2;
                    half_max_x = eP_profs.half_max_x;
                    t_prof_fit = out_t_prof.yFit_ups;
                    t_prof_fit_txt = ...
                        sprintf('fit with %s', out_t_prof.fitModel);
                    half_max_t_1 = eP_profs.half_max_t_1;
                    half_max_t_2 = eP_profs.half_max_t_2;
                    half_max_t = eP_profs.half_max_t;
                    t_max = eP_profs.t_max;
                    bs_t = eP_profs.bs_t;
                    v_max = eP_profs.v_max;
                    t0 = eP_profs.t0;
                case 'estimate from event img'
                    txt_a = ['parameters from profiles', ...
                        '(1 um and 5 ms wide) crossing peak of event',...
                        '(from image smoothed with 2D spline)'];
                    x_ups = x;
                    x_prof_fit = x_event_profS;
                    x_prof_fit_txt = ...
                        'profile from image smoothed with 2D spline';
                    if ~exist('eParamsEst','var')
                        switch calcMethodPrev
                            case 'peakXTProfiles'
                                eParamsEst = eP_profs;
                            case '2DGauss'
                                eParamsEst = eP_2Dfit;
                        end
                    end
                    half_max_x_1 = eParamsEst.half_max_x_1;
                    half_max_x_2 = eParamsEst.half_max_x_2;
                    half_max_x = eParamsEst.half_max_x;
                    t_ups = t;
                    t_prof_fit = t_event_profS;
                    t_prof_fit_txt = ...
                        'profile from image smoothed with 2D spline';
                    half_max_t_1 = eParamsEst.half_max_t_1;
                    half_max_t_2 = eParamsEst.half_max_t_2;
                    half_max_t = eParamsEst.half_max_t;
                    t_max = eParamsEst.t_max;
                    try 
                        bs_t = eParamsEst.bs_t;
                    catch
                        bs_t = eParamsEst.bs;
                    end
                    v_max = eParamsEst.v_max;
                    t0 = eParamsEst.t0;    
            end

            % 1.text axes
            ax_txt_a = axes(hf, 'Units','normalized', ...
                'Position',[dx 1-dy 1-2*dx dy], ...
                'Visible','off');
            text(ax_txt_a, 0, 0.5,  txt_a, 'FontUnits','normalized', ...
                'FontSize',0.4, 'FontWeight','bold', ...
                'HorizontalAlignment','left')
            % image axes
            ha1 = axes(hf, 'Units','normalized', ...
                'Position',[dx 1-h_a-2*dy w_a*1.5 h_a]);
            switch calcMethod
                case '2DGauss'
                    mesh(T, X, imgE, ...
                        'LineStyle','-', 'LineWidth',0.5, ...
                        'FaceColor','none', 'EdgeColor','k', ...
                        'EdgeAlpha',0.4, 'Parent',ha1)
                    %plot3(T, X, D, 'LineStyle','-', 
                    % 'Marker','.', 'Color',[0.1 0.1 0.1], 'Parent',ha1)
                    hold on
                    surf(T, X, imgE_fit, ...
                        'FaceAlpha',0.6, 'EdgeColor','none', ...
                        'FaceColor','interp', 'Parent',ha1)
                    line(t_ups, X_ups(rFitMaxUps,:), t_prof_fit, ...
                        'Parent',ha1, 'LineWidth',3, 'Color','k')
                    line(T_ups(:,cFitMaxUps), x_ups, x_prof_fit, ...
                        'Parent',ha1, 'LineWidth',3, 'Color','k')
                    set(ha1, 'XLim',[min(t) max(t)], ...
                        'YLim',[min(x) max(x)], ...
                        'ZLim',getAxisLimits(imgE, 5))
                    colormap(jet)
                    view(-15,40)
                    set(ha1,'YDir','normal', 'FontSize',fontSzNum)
                    title(ha1, [fig_title_txt,' -- 2D gauss fit'], ...
                        'FontSize',fontSzT)
                    xlabel(ha1,'t (ms)', 'FontSize',fontSzL)
                    ylabel(ha1,'x (\mum)', 'FontSize',fontSzL)
                    zlabel(hObjs.ax_prof.YLabel.String, 'FontSize',fontSzL)
                otherwise
                    image(imgE, 'YData',[min(x_ups) max(x_ups)], ...
                        'XData',[min(t_ups) max(t_ups)], ...
                        'CDataMapping','scaled', 'Parent',ha1)
                    % show lines of areas from where profiles are
                    % calculated
                    line(get(gca,'XLim'), ...
                        [(r_m-1-(n_px_t-1)/2)*pxSzX (r_m-1-(n_px_t-1)/2)*pxSzX], ...
                        'Parent',ha1, 'LineWidth',1, ...
                        'Color','k', 'LineStyle','-');
                    line(get(gca,'XLim'), ...
                        [(r_m-1+(n_px_t-1)/2)*pxSzX (r_m-1+(n_px_t-1)/2)*pxSzX], ...
                        'Parent',ha1, 'LineWidth',1, ...
                        'Color','k', 'LineStyle','-');
                    line([(c_m-1-(n_px_x-1)/2)*pxSzT (c_m-1-(n_px_x-1)/2)*pxSzT], ...
                        get(gca,'YLim'), 'Parent',ha1, ...
                        'LineWidth',1, 'Color','k', 'LineStyle','-');
                    line([(c_m-1+(n_px_x-1)/2)*pxSzT (c_m-1+(n_px_x-1)/2)*pxSzT], ...
                        get(gca,'YLim'), 'Parent',ha1, ...
                        'LineWidth',1, 'Color','k', 'LineStyle','-');
                    set(ha1, 'FontSize',fontSzNum)
                    old_img_cLims = get(ha1,'Clim');
                    % in case taht in image of events there is also a part
                    % of another event with much higher amplitude
                    set(ha1, 'Clim', [old_img_cLims(1) ...
                        prctile(imgE(imgE_m),99,'all')])
                    title(ha1, fig_title_txt, 'FontSize',fontSzT)
                    xlabel(ha1,'t (ms)', 'FontSize',fontSzL)
                    ylabel(ha1,'x (\mum)', 'FontSize',fontSzL)
            end

            % x-profile axes
            ha2 = axes(hf, 'Units','normalized', ...
                'Position',[ha1.Position(1)+ha1.Position(3)+dx ...
                            1-h_a-2*dy ...
                            w_a h_a]);
            line(x, x_event_prof, 'Parent',ha2, ...
                'LineWidth',1, 'Color','k', 'Display','data')
            line(x_ups, x_prof_fit, 'Parent',ha2, ...
                'LineWidth',2, 'Color','r', 'Display',x_prof_fit_txt)
            set(ha2,'FontSize',fontSzNum)
            title(ha2, ax_profX_title_txt, 'FontSize',fontSzT)
            ylabel(ha2, hObjs.ax_prof.YLabel.String, 'FontSize',fontSzL)
            xlabel(ha2, 'x (\mum)', 'FontSize',fontSzL)
            xlim(ha2, [min(x_ups) max(x_ups)])
            ylim(ha2, getAxisLimits([x_event_prof(:);x_prof_fit(:)], 5))
            line([half_max_x_1 half_max_x_2], ...
                 [half_max_x half_max_x], 'Parent',ha2,...
                'LineWidth',2, 'Color','b', 'Display','FWHM')
            hl_2 = legend(ha2,'show');
            hl_2.Location = 'best';
            hl_2.FontSize = fontSzLegend;
           
            % t-profile axes
            ha3 = axes(hf, 'Units','normalized', ...
                'Position',[ha2.Position(1)+ha2.Position(3)+dx ...
                            1-h_a-2*dy ...
                            w_a h_a]);
            line(t, t_event_prof, 'Parent',ha3, ...
                'LineWidth',1, 'Color','k', 'Display','data')
            line(t_ups, t_prof_fit, 'Parent',ha3, ...
                'LineWidth',2, 'Color','r', 'Display',t_prof_fit_txt)
            set(ha3,'FontSize',fontSzNum)
            title(ha3, ax_profT_title_txt, 'FontSize',fontSzT)
            xlabel(ha3, 't (ms)', 'FontSize',fontSzL)
            ylabel(ha3, hObjs.ax_prof.YLabel.String, 'FontSize',fontSzL)
            xlim(ha3, [min(t_ups) max(t_ups)])
            ylim(ha3, getAxisLimits([t_event_prof(:);t_prof_fit(:)], 5))
            line([half_max_t_1 half_max_t_2], ...
                 [half_max_t half_max_t], 'Parent',ha3,...
                'LineWidth',2, 'Color','b', 'Display','FDHM')
            line([t_max t_max], [bs_t v_max], 'Parent',ha3, ...
                'LineWidth',2, 'Color','g', 'Display','amplitude')
            line([t0 t_max], [bs_t bs_t], 'Parent',ha3, ...
                'LineWidth',2, 'Color','m', 'Display','TTP')
            hl_3 = legend(ha3, 'show');
            hl_3.Location = 'best';
            hl_3.FontSize = fontSzLegend;

            % 2.text axes, parameters of analyzed event
            params_txt = {
                [sprintf('amplitude = %0.2f',eventParams.amplitude(i,1)), ...
                ' (',char(916),'F/F_0)'],...
                sprintf('TTP = %0.2f (ms)',eventParams.TTP(i,1)), ...
                sprintf('FDHM = %0.2f (ms)',eventParams.FDHM(i,1)), ...
                [sprintf('FWHM = %0.2f ',eventParams.FWHM(i,1)), ...
                '(',char(181),'m)'],...
                [sprintf('sparkMass = %0.2f ',eventParams.sparkMass(i,1)), ...
                '(',char(916),'F/F_0*',char(181),'m^3)']};
            ax_txt_b = axes(hf, 'Units','normalized', ...
                'Position',[dx 1-4*dy-h_a-2*dy 1-2*dx 2*dy], ...
                'Visible','off');
            text(ax_txt_b, 0, 1,  ...
                [txt_b, ...
                join(params_txt, ' // ')], ...
                'FontUnits','pixels', ...
                'FontSize',18, 'FontWeight','normal', ...
                'HorizontalAlignment','left', ...
                'VerticalAlignment','top')               
        end
        
        % if method for parameters estimation was changed,
        % due to error in selected method, 
        % then after using simpler method for current iteration calculation
        % change it back to selected one
        calcMethod = calcMethodPrev;
        
        clearvars D t x t_prof x_prof t_ups x_ups ...
            t_prof_ups x_prof_ups y_t1 y_t2 D1 y_x1 y_x2
        
        clearvars rows_e cols_e eP_profs eP_2Dfit ...
            out_t_prof out_x_prof p0 t_spark_prof x_spark_prof
        
        % Report current estimate in the waitbar's message field
        waitbar(i/length(statSparks), hw, ...
            sprintf('%d %%',round(100*i/length(statSparks))))
        
    end
    
    % delete waitbar
    delete(hw), clearvars hw
    
else
    eventParams.amplitude = [];
    eventParams.t0 = [];
    eventParams.TTP = [];
    eventParams.FDHM = [];
    eventParams.FWHM = [];
    eventParams.sparkMass = [];
    eventParams.tauD = [];
    eventParams.AUC_2DFit = [];
    eventParams.sparkROINum = [];
end

% put back pointer to arrow
set(mainFig,'Pointer','arrow')
drawnow

end



