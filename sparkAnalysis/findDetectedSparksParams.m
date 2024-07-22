function [eventParams, analyzedEvntsBrowserTbl] = findDetectedSparksParams(...
    img, statSparks, mainFig ,calcMethod, sparksRec, posOfROI, ...
    startOfSpark, endOfSpark, prevFitCoeffs, normImgFlag, ...
    bsDetSensitivity, smoothSpan)

% img = filtered data or filtered normalized data
% statSparks = statistic from sparks detection
% startOfSpark, endOfSpark, prevFitCoeffs = previous fit with fitEventRise
% function

% posOfROI = from spark recovery analysis
if isempty(posOfROI)
    posOfROI = nan;
end

% get data
imgData = getappdata(mainFig,'imgData');
pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;
hObjs = getappdata(mainFig,'hObjs');

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
    % table to show analyzed events
    analyzedEvntsBrowserTbl = table();
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
            prevSparkFitCoeffs, ...
            bsDetSensitivity, smoothSpan);
        
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
                    % save parameters
                    eP_2Dfit_fiedN = fieldnames(eP_2Dfit);
                    eventParams.calcMethod{i,1} = calcMethod;
                    eventParams.amplitude(i,1) = eP_2Dfit.amplitude;
                    eventParams.t0(i,1) = eP_2Dfit.t0 + eventROIstart_t;
                    eventParams.(eP_2Dfit_fiedN{contains(eP_2Dfit_fiedN, 't0_lineFit_')})(i,1) = ...
                        eP_2Dfit.(eP_2Dfit_fiedN{contains(eP_2Dfit_fiedN, 't0_lineFit_')}) ...
                        + eventROIstart_t;
                    eventParams.(eP_2Dfit_fiedN{contains(eP_2Dfit_fiedN, 'TTP_lineFit_')})(i,1) = ...
                        eP_2Dfit.(eP_2Dfit_fiedN{contains(eP_2Dfit_fiedN, 'TTP_lineFit_')});
                    eventParams.TTP(i,1) = eP_2Dfit.TTP;
                    eventParams.FDHM(i,1) = eP_2Dfit.FDHM;
                    eventParams.FWHM(i,1) = eP_2Dfit.FWHM;
                    eventParams.tauD(i,1) = nan;
                    eventParams.tauR(i,1) = nan;
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
                    %keyboard
                    eventParams.calcMethod{i,1} = calcMethod;
                    eventParams.amplitude(i,1) = nan;
                    eventParams.t0(i,1) = nan;
                    eventParams.t0_fittedLine(i,1) = nan;
                    eventParams.TTP_fittedLine(i,1) = nan;
                    eventParams.TTP(i,1) = nan;
                    eventParams.FDHM(i,1) = nan;
                    eventParams.FWHM(i,1) = nan;
                    eventParams.sparkMass(i,1) = nan;
                    eventParams.tauD(i,1) = nan;
                    eventParams.tauR(i,1) = nan;
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
                        out_t_prof = fitOneWholeEvent(t, t_event_prof, ...
                            "constThenSpline", t0=t0, bs=bs, ...
                            x_ups=t_ups, mProf=tProf_m);
                        bs_t_prof = out_t_prof.bs;
                        t0_t_prof = out_t_prof.t0;
                        tauD = out_t_prof.tauD;
                        tauR = out_t_prof.tauR;
                        % fit gaussian or exp modified gaussian or spline
                        % to spatial profile
                        out_x_prof = fitOneWholeEvent(x, x_event_prof, ...
                            "Gauss", bs=bs, ...
                            x_ups=x_ups, mProf=xProf_m, ...
                            addMoreWeightToEvent=true);
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
                        if isempty(loc_ups)
                            [~, loc_ups] = max(out_t_prof.yFit_ups(tProf_m_ups));
                            loc_ups = t_ups(loc_ups + sE_ups);
                        end
                        [~, peakData_ups.pos] = ...
                            min( abs( t_ups - loc_ups ) );
                        peakData_ups.val = ...
                            out_t_prof.yFit_ups(peakData_ups.pos);
                        % get parameters from profile fits
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
                        % estimate tauD and tauR
                        out_t_prof = fitOneWholeEvent(t, t_event_profS, ...
                            "spline", t0=eP_profs.t0, bs=eP_profs.bs, ...
                            mProf=tProf_m);
                        tauD = out_t_prof.tauD;
                        tauR = out_t_prof.tauR;
                    end
                    % save parameters
                    eP_profs_fiedN = fieldnames(eP_profs);

                    eventParams.calcMethod{i,1} = calcMethod;
                    eventParams.amplitude(i,1) = eP_profs.amplitude;
                    eventParams.t0(i,1) = eP_profs.t0 + eventROIstart_t;
                    eventParams.(eP_profs_fiedN{contains(eP_profs_fiedN, 't0_lineFit_')})(i,1) = ...
                        eP_profs.(eP_profs_fiedN{contains(eP_profs_fiedN, 't0_lineFit_')}) ...
                        + eventROIstart_t;
                    eventParams.(eP_profs_fiedN{contains(eP_profs_fiedN, 'TTP_lineFit_')})(i,1) = ...
                        eP_profs.(eP_profs_fiedN{contains(eP_profs_fiedN, 'TTP_lineFit_')});
                    eventParams.TTP(i,1) = eP_profs.TTP;
                    eventParams.FDHM(i,1) = eP_profs.FDHM;
                    eventParams.FWHM(i,1) = eP_profs.FWHM;
                    eventParams.sparkMass(i,1) = eP_profs.sparkMass;
                    eventParams.tauD(i,1) = tauD;
                    eventParams.tauR(i,1) = tauR;
                    eventParams.AUC_2DFit(i,1) = nan;
                catch
                    %keyboard
                    % if cannot get profiles for some reason
                    eventParams.calcMethod{i,1} = calcMethod;
                    eventParams.amplitude(i,1) = nan;
                    eventParams.t0(i,1) = nan;
                    eventParams.t0_fittedLine(i,1) = nan;
                    eventParams.TTP_fittedLine(i,1) = nan;
                    eventParams.TTP(i,1) = nan;
                    eventParams.FDHM(i,1) = nan;
                    eventParams.FWHM(i,1) = nan;
                    eventParams.sparkMass(i,1) = nan;
                    eventParams.tauD(i,1) = nan;
                    eventParams.tauR(i,1) = nan;
                    eventParams.AUC_2DFit(i,1) = nan;
                end
        
                
            case 'estimate from event img'
                try
                    % get parameters of spark from profiles without 2D fitting
                    % use coeficient from previous fitting (t0, tauR, A, baseline)
                    % get params from image of event
                    if ~isempty(t0)
                        eParamsEst = getParametersOfEventProfile( ...
                            t, t_event_profS, ...
                            x, x_event_profS, ...
                            bs, bs, t0, ...
                            imgData.blank, ...
                            peakData_tProf, tProf_m, normImgFlag);
                        % estimate tauD and tauR
                        out_t_prof = fitOneWholeEvent(t, t_event_profS, ...
                            "spline", t0=t0, bs=bs, ...
                            mProf=tProf_m);
                        tauD = out_t_prof.tauD;
                        tauR = out_t_prof.tauR;
                    else
                        eParamsEst = ...
                            getParametersOfEventFromRawImage(...
                            imgE, imgEs, imgE_m, ...
                            [], [], r_m, c_m, ...
                            t_event_profS, x_event_profS, t, x, ...
                            imgData, normImgFlag);
                        % estimate tauD and tauR
                        out_t_prof = fitOneWholeEvent(t, t_event_profS, ...
                            "spline", t0=eParamsEst.t0, bs=eParamsEst.bs, ...
                            mProf=tProf_m);
                        tauD = out_t_prof.tauD;
                        tauR = out_t_prof.tauR;
                    end
                    eParamsEst_fiedN = fieldnames(eParamsEst);

                    eventParams.calcMethod{i,1} = calcMethod;
                    eventParams.amplitude(i,1) = eParamsEst.amplitude;
                    eventParams.t0(i,1) = eParamsEst.t0 + eventROIstart_t;
                    eventParams.(eParamsEst_fiedN{contains(eParamsEst_fiedN, 't0_lineFit_')})(i,1) = ...
                        eParamsEst.(eParamsEst_fiedN{contains(eParamsEst_fiedN, 't0_lineFit_')}) ...
                        + eventROIstart_t;
                    eventParams.(eParamsEst_fiedN{contains(eParamsEst_fiedN, 'TTP_lineFit_')})(i,1) = ...
                        eParamsEst.(eParamsEst_fiedN{contains(eParamsEst_fiedN, 'TTP_lineFit_')});
                    eventParams.TTP(i,1) = eParamsEst.TTP;
                    eventParams.FDHM(i,1) = eParamsEst.FDHM;
                    eventParams.FWHM(i,1) = eParamsEst.FWHM;
                    eventParams.sparkMass(i,1) = eParamsEst.sparkMass;
                    eventParams.tauD(i,1) = tauD;
                    eventParams.tauR(i,1) = tauR;
                    eventParams.AUC_2DFit(i,1) = nan;
                catch
                    % if cannot get estimate of params for some reason
                    %keyboard
                    eventParams.calcMethod{i,1} = calcMethod;
                    eventParams.amplitude(i,1) = nan;
                    eventParams.t0(i,1) = nan;
                    eventParams.t0_fittedLine(i,1) = nan;
                    eventParams.TTP_fittedLine(i,1) = nan;
                    eventParams.TTP(i,1) = nan;
                    eventParams.FDHM(i,1) = nan;
                    eventParams.FWHM(i,1) = nan;
                    eventParams.sparkMass(i,1) = nan;
                    eventParams.tauD(i,1) = nan;
                    eventParams.tauR(i,1) = nan;
                    eventParams.AUC_2DFit(i,1) = nan;
                end  
        end
        
        eventParams.sparkROINum(i,1) = sparkROInum;
        eventParams.BoundingBox_x(i,1) = statSparks(i).BoundingBox(1);
        eventParams.BoundingBox_y(i,1) = statSparks(i).BoundingBox(2);
        eventParams.BoundingBox_w(i,1) = statSparks(i).BoundingBox(3);
        eventParams.BoundingBox_h(i,1) = statSparks(i).BoundingBox(4);
        
        % create table with all data neccesary to show analyzed events
        % setup titles of figures and axes
        switch getappdata(mainFig, 'analysisType')
            case 'spark detection'
                fig_title_txt = ...
                    sprintf('spark from ROI #: %d',sparkROInum);
                ax_profX_title_txt = 'x profile';
                ax_profT_title_txt = 't profile';
                txt_b = {sprintf('parameters of spark from ROI #%d', ...
                    sparkROInum)};
            case 'spark recovery ryanodine'
                fig_title_txt = ...
                    sprintf('posOfRoi %.2f \\mum; spark #%d', ...
                    posOfROI*pxSzX, sparkROInum);
                ax_profX_title_txt = 'x profile';
                ax_profT_title_txt = 't profile';
                    % sprintf('posOfRoi %.2f \\mum; spark #%d -- t profile', ...
                    % posOfROI*pxSzX, sparkROInum);
                txt_b = {[sprintf('posOfRoi %.2f ',posOfROI*pxSzX), ...
                    char(181),'m -- ', sprintf('parameters of spark #%d', ...
                    sparkROInum)]};
        end
        % save data to struct
        axDesc = struct('fig_title_txt',fig_title_txt,...
            'ax_profX_title_txt',ax_profX_title_txt, ...
            'ax_profT_title_txt',ax_profT_title_txt, ...
            'txt_b',txt_b, ...
            'ylabel',hObjs.ax_prof.YLabel.String);
        % data to show estimated parameters and profiles
        switch calcMethod
            case '2DGauss'
                txt_a = ['2D exponentially modified gaussian fit ', ...
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
                % show fits for exp rise and decay
                t_rise = [];
                t_profR_fit = [];
                t_decay = [];
                t_profD_fit = [];
                T = T;
                X = X;
                x_ups_peak = X_ups(rFitMaxUps,:);
                t_ups_peak = T_ups(:,cFitMaxUps);
                imgE_fit = imgE_fit;

            case 'peakXTProfiles'
                txt_a = '1 µm and 5 ms wide peak profiles';
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
                % show fits for exp rise and decay
                if exist("peakData_ups", "var")
                    t_rise = t_ups(1:peakData_ups.pos);
                    t_decay = t_ups(peakData_ups.pos:end);
                else
                    t_rise = t(1:out_t_prof.evntMaxPosPx);
                    t_decay = t(out_t_prof.evntMaxPosPx:end);
                end
                if ~isempty(out_t_prof.coeffRise)
                    t_profR_fit = out_t_prof.yFitFunR( ...
                        out_t_prof.coeffRise, t_rise);
                else
                    t_profR_fit = nan(size(t_rise));
                end
                if ~isempty(out_t_prof.coeffDecay)
                    t_profD_fit = out_t_prof.yFitFunD( ...
                        out_t_prof.coeffDecay, t_decay);
                else
                    t_profD_fit = nan(size(t_decay));
                end
                T = [];
                X = [];
                x_ups_peak = [];
                t_ups_peak = [];
                imgE_fit = [];

            case 'estimate from event img'
                txt_a = ['1 µm and 5 ms wide peak profiles ',...
                    '(image smoothed with 2D spline)'];
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
                % show fits for exp rise and decay
                t_rise = t(1:out_t_prof.evntMaxPosPx);
                t_profR_fit = out_t_prof.yFitFunR( ...
                    out_t_prof.coeffRise, t_rise);
                t_decay = t(out_t_prof.evntMaxPosPx:end);
                t_profD_fit = out_t_prof.yFitFunD( ...
                    out_t_prof.coeffDecay, t_decay);
                T = [];
                X = [];
                x_ups_peak = [];
                t_ups_peak = [];
                imgE_fit = [];
        end
        % save data to struct
        axData = struct('txt_a',txt_a,...
            'imgE', imgE, ...
            'imgE_m', imgE_m, ...
            'x', x, ...
            'x_ups',x_ups, ...
            'x_event_prof',x_event_prof, ...
            'x_prof_fit',x_prof_fit, ...
            'x_prof_fit_txt',x_prof_fit_txt, ...
            'half_max_x_1',half_max_x_1, ...
            'half_max_x_2',half_max_x_2, ...
            'half_max_x',half_max_x, ...
            't', t, ...
            't_ups',t_ups, ...
            't_event_prof', t_event_prof, ...
            't_prof_fit',t_prof_fit, ...
            't_prof_fit_txt',t_prof_fit_txt, ...
            'half_max_t_1',half_max_t_1, ...
            'half_max_t_2',half_max_t_2, ...
            'half_max_t',half_max_t, ...
            't_max',t_max, ...
            'bs_t',bs_t, ...
            'v_max',v_max, ...
            't0',t0, ...
            't_rise',t_rise, ...
            't_profR_fit',t_profR_fit, ...
            't_decay',t_decay, ...
            't_profD_fit',t_profD_fit, ...
            'r_m', r_m, ...
            'n_px_t', n_px_t, ...
            'c_m', c_m, ...
            'n_px_x', n_px_x, ...
            'T', T, ...
            'X', X, ...
            'x_ups_peak', x_ups_peak, ...
            't_ups_peak', t_ups_peak, ...
            'imgE_fit', imgE_fit);

        % save data to table
        analyzedEvntsBrowserTbl = [ ...
            analyzedEvntsBrowserTbl; 
            table(sparkROInum, {getappdata(mainFig, 'analysisType')}, ...
                {calcMethod}, {axDesc}, {axData}, ...
                {structfun(@(x) x(i),eventParams, "UniformOutput",false)}, ...
                'VariableNames',{'sparkROInum', 'analysisType', ...
                'calcMethod', 'axDesc', 'axData', 'evntParams'})
                ];

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
    eventParams.calcMethod = [];
    eventParams.amplitude = [];
    eventParams.t0 = [];
    eventParams.t0_lineFit_10_90 = [];
    eventParams.TTP = [];
    eventParams.TTP_lineFit_10_90 = [];
    eventParams.FDHM = [];
    eventParams.FWHM = [];
    eventParams.sparkMass = [];
    eventParams.tauD = [];
    eventParams.tauR = [];
    eventParams.AUC_2DFit = [];
    eventParams.sparkROINum = [];
    eventParams.BoundingBox_x = [];
    eventParams.BoundingBox_y = [];
    eventParams.BoundingBox_w = [];
    eventParams.BoundingBox_h = [];

    analyzedEvntsBrowserTbl = [];
end

% put back pointer to arrow
set(mainFig,'Pointer','arrow')
drawnow

end



