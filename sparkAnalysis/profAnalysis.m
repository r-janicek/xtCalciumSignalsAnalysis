function profAnalysis(~,~,mainFig)
% analyze all selected profiles
% get data 
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');
profileAnalysis = getappdata(mainFig,'profileAnalysis');

if ~isfield(profileAnalysis,'selectedROIs')
   return 
end

h = findobj('-regexp','Name','repetitive sparks analysis');
if ~isempty(h)
    close(h)
end 

% mouse cursor
set(mainFig,'Pointer','watch')
drawnow

selectedROIs = profileAnalysis.selectedROIs;
t = imgData.t;
pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;

fit_result = num2cell(zeros(height(selectedROIs),2));
eventsParameters = num2cell(zeros(height(selectedROIs),1));

for i=1:height(selectedROIs)
    
    if any(strcmp(selectedROIs.Properties.VariableNames,'wholeProfileFit'))
        % fitted profile
        wholeProfileFit = selectedROIs.wholeProfileFit{i,1};
        % normalized profile
        prof = wholeProfileFit.yN;
        
        locs = cell2mat(wholeProfileFit.eventsPeaks(:,2));
        
        pks = prof(round(locs./pxSzT)+1);
        
        coefPrevFit = selectedROIs.prevFitCoef{i,1};   
        coefPrevFit(:,3) = pks;
        coefPrevFit(:,4) = coefPrevFit(:,4)./ ...
            mean(wholeProfileFit.y(wholeProfileFit.baselineM)) - 1;   
    else
        % normalized profile
        prof = mean(selectedROIs.dataProfile{i,1},1);
        coefPrevFit = selectedROIs.prevFitCoef{i,1};
        locs = cell2mat(selectedROIs.eventsPeaks{i,1}(:,2));    
        pks = prof(round(locs./pxSzT)+1);
    end
          
    dataWholeArea = selectedROIs.dataWholeArea{i,1};  
    dataWholeAreaRaw = selectedROIs.dataWholeAreaRaw{i,1}; 
    statEvents = selectedROIs.regionsProp{i,1};   
    posOfROI = selectedROIs.positionOfRoi(i,1);
    smooth_span = selectedROIs.smooth_span(i,1);
    bs_crit = selectedROIs.bs_crit(i,1);
    sSp = selectedROIs.startPosOfSpark{i,1};
    eSp = selectedROIs.endPosOfSpark{i,1};
    prof_t_evnts_m = false(numel(prof),1);
    % create events mask for profile
    for e=1:length(statEvents)
        prof_t_evnts_m(statEvents(e).SubarrayIdx{2}, 1) = true;
    end
    % fit only rise of spark fun(t0,u,tR,A,y0)
    evntRiseFits = fitEventRise(t, prof, "local", ...
        peaks_vals=pks, peaks_locs=locs, ...
        coefPrevFit=coefPrevFit, ...
        sSpPrev=sSp, eSpPrev=eSp, fitTol=1e-9, ...
        bs_crit=bs_crit, smooth_span=smooth_span, ...
        evntsMask=prof_t_evnts_m, numOfFitIter=1000);

    fit_result(i,1) = {evntRiseFits.coef};
    fit_result(i,2) = {evntRiseFits.sp_fit};
     
    % calculate sparks parameters
    if hObjs.check_sparkParams.Value
        % selected method of parameters calculation
        switch hObjs.calcSpParamMethodRBgroup.SelectedObject.String
            case '2D gaussian fit'
                calcMethod = '2DGauss'; 
            case 'max crossing profiles'
                calcMethod = 'peakXTProfiles';
        end
        % check if use normalized image or filtered raw image to calculate spark
        % parameters
        if hObjs.check_useNormalizedImg.Value
            img = dataWholeArea;
        else
            % filter raw image
            img = imgFiltering(dataWholeAreaRaw, pxSzT, pxSzX);
        end
        % calculate parameters of selected events
        eventsParams_roi = findDetectedSparksParams(img, ...
            statEvents, mainFig, calcMethod, [], posOfROI, ...
            evntRiseFits.startOfEvent, evntRiseFits.endOfEvent, ...
            evntRiseFits.coef, hObjs.check_useNormalizedImg.Value);
        eventsParameters(i,1) = {eventsParams_roi};

    else
        eventsParams_roi.amplitude = nan(size(statEvents));
        eventsParams_roi.TTP = nan(size(statEvents));
        eventsParams_roi.TTP_fittedLine = nan(size(statEvents));
        eventsParams_roi.FDHM = nan(size(statEvents));
        eventsParams_roi.FWHM = nan(size(statEvents));
        eventsParams_roi.sparkMass = nan(size(statEvents));
        eventsParams_roi.tauD = nan(size(statEvents));
        eventsParams_roi.AUC_2DFit = nan(size(statEvents));
        
        eventsParameters(i,1) = {eventsParams_roi}; 
    end
end

% create result structure
selectedROIs.finalFitResults = fit_result;
selectedROIs.eventsParameters = eventsParameters;

% analyze peaks ratio and distances
for i = 1:height(selectedROIs)
    
    if ~isempty(selectedROIs.eventsPeaks{i,1})
        
        % calculate parameters of events from whole profile fit
        if any(strcmp(selectedROIs.Properties.VariableNames,'wholeProfileFit'))
            paramsFittedEvents = ...
                calcParametersOfEventsFromWholeProfileFit(...
                    t, selectedROIs.wholeProfileFit{i,1}.profFit);
            % normalized profiles
            profN = selectedROIs.wholeProfileFit{i,1}.yN;
            profN = profN(:);
            
        else            
            names_params = {'t0 (ms)','tPeak (ms)','Apeak (F/F0)','TTP (ms)','FDHM (ms)','AUC (ms*F/F0)','tauD (ms)','firstDerMaxVal (deltaF/F0)','firstDerMaxValPos (ms)'};
            paramsFittedEvents = [names_params;num2cell(nan(size(selectedROIs.eventsPeaks{i,1},1),9))];
            % normalized profiles
            profN = mean(selectedROIs.dataProfile{i,1},1);
            profN = profN(:);
        end
        
        selectedROIs.normProf(i) = {profN};
        
        %detected peak amplitude minus baseline (y0 from fitting of rise of spark)
        y0 = selectedROIs.finalFitResults{i,1}(:,4);
        peaks_P = cell2mat(selectedROIs.eventsPeaks{i,1}(:,2));
        peaks_A = profN(round(peaks_P./pxSzT)+1)-y0;
        
        
        % calculate 1. derivative max amplitude and position for
        % all selected events in normalized profile
        [firstDerMaxValPos,firstDerMaxVal] = getMaxOfFirstDerOfEvents(profN, selectedROIs.startPosOfSpark{i}, round(peaks_P./pxSzT)+1);
        firstDerMaxValPos = t(firstDerMaxValPos); % in ms
        firstDerMaxValPos = firstDerMaxValPos(:);
        
        % save peaks values
        selectedROIs.eventsPeaks{i,1}(:,1) = num2cell(peaks_A);
        
        sparkParams = [selectedROIs.eventsParameters{i,1}];
        sparkAmp = sparkParams.amplitude;
        sparkMass = sparkParams.sparkMass;
       
        peaksWholeProfFit_A = [paramsFittedEvents{2:end,3}]';
        peaksWholeProfFit_P = [paramsFittedEvents{2:end,2}]';
        peaksWholeProfFit_TTP = [paramsFittedEvents{2:end,4}]';
        peaksWholeProfFit_FDHM = [paramsFittedEvents{2:end,5}]';
        peaksWholeProfFit_AUC = [paramsFittedEvents{2:end,6}]';
        peaksWholeProfFit_tauD = [paramsFittedEvents{2:end,7}]';
        peaksWholeProfFit_derA = [paramsFittedEvents{2:end,8}]'; % 1. derivative
        peaksWholeProfFit_derP = [paramsFittedEvents{2:end,9}]';
        

        if numel(peaks_P) > 1
            % peaks positions difference
            peaks_P_diff = diff(peaks_P);

            peaks_derP_diff = diff(firstDerMaxValPos);
            peaks_derP_diff = peaks_derP_diff(:);

            peaksWholeProfFit_derP_diff = diff(peaksWholeProfFit_derP);


            for j = 1:length(peaks_P_diff)

                % ratios/diff of data from profile or spark
                peaks_A_ratio(j,1) = peaks_A(j+1) / peaks_A(j);

                peaks_derA_ratio(j,1) = firstDerMaxVal(j+1) / firstDerMaxVal(j);

                sparkMassRatio(j,1) = sparkMass(j+1) / sparkMass(j);
                sparkAmpRatio(j,1) = sparkAmp(j+1) / sparkAmp(j);

                % ratios/diff of data from fit of profile
                peaksWholeProfFit_A_ratio(j,1) = peaksWholeProfFit_A(j+1) / peaksWholeProfFit_A(j);
                peaksWholeProfFit_P_diff(j,1) = peaksWholeProfFit_P(j+1) - peaksWholeProfFit_P(j);
                peaksWholeProfFit_TTP_ratio(j,1) = peaksWholeProfFit_TTP(j+1) / peaksWholeProfFit_TTP(j);
                peaksWholeProfFit_FDHM_ratio(j,1) = peaksWholeProfFit_FDHM(j+1) / peaksWholeProfFit_FDHM(j);
                peaksWholeProfFit_AUC_ratio(j,1) = peaksWholeProfFit_AUC(j+1) / peaksWholeProfFit_AUC(j);
                peaksWholeProfFit_tauD_ratio(j,1) = peaksWholeProfFit_tauD(j+1) / peaksWholeProfFit_tauD(j);

                peaksWholeProfFit_derA_ratio(j,1) = peaksWholeProfFit_derA(j+1) / peaksWholeProfFit_derA(j);

                % decide whether accept pair or not
                if j==1

                    if peaks_P(1) < 200 % 200 ms taken from Sobie analysis
                        accepted(j,1) = false;
                    else
                        accepted(j,1) = true;
                    end

                else

                    if (peaks_P_diff(j-1) < 200)

                        accepted(j,1) = false;

                    else
                        accepted(j,1) = true;

                    end

                end

            end

        else
            peaks_A_ratio = [];
            peaks_P_diff = [];
            accepted = [];
            sparkMassRatio = [];
            sparkAmpRatio = [];
            
            peaks_derP_diff = [];
            peaks_derA_ratio = [];
            
            peaksWholeProfFit_A_ratio = [];
            peaksWholeProfFit_P_diff = [];
            peaksWholeProfFit_TTP_ratio = [];
            peaksWholeProfFit_FDHM_ratio = [];
            peaksWholeProfFit_AUC_ratio = [];
            peaksWholeProfFit_tauD_ratio = [];
            
            peaksWholeProfFit_derP_diff = [];
            peaksWholeProfFit_derA_ratio = [];
            
        end
        
    else
        peaks_A_ratio = [];
        peaks_P_diff = [];
        accepted = [];
        sparkMassRatio = [];
        sparkAmpRatio = [];
        
        peaks_derP_diff = [];
        peaks_derA_ratio = [];
             
        peaksWholeProfFit_A_ratio = [];
        peaksWholeProfFit_P_diff = [];
        peaksWholeProfFit_TTP_ratio = [];
        peaksWholeProfFit_FDHM_ratio = [];
        peaksWholeProfFit_AUC_ratio = [];
        peaksWholeProfFit_tauD_ratio = [];
        
        peaksWholeProfFit_derP_diff = [];
        peaksWholeProfFit_derA_ratio = [];
        
    end
   
    % save results
    T_result = table(peaks_A_ratio,peaks_P_diff,accepted,sparkAmpRatio,sparkMassRatio,...
        peaksWholeProfFit_A_ratio,peaksWholeProfFit_P_diff,peaksWholeProfFit_TTP_ratio,...
        peaksWholeProfFit_FDHM_ratio,peaksWholeProfFit_AUC_ratio,peaksWholeProfFit_tauD_ratio,...
        peaks_derP_diff,peaks_derA_ratio,peaksWholeProfFit_derP_diff,peaksWholeProfFit_derA_ratio,...
        'VariableNames',{'AmplitudeRatio' 'peakPosDiff' 'acceptedPair' 'sparkAmpRatio' 'sparkMassRatio' ...
        'AmplRatioFit' 'peakPosDiffFit' 'TTPratioFit' 'FDHMratioFit' 'AUCratio' 'tauDratioFit' ...
        'peakPosDiffDeriv','amplDerivRatio','peakPosDiffDerivFit','amplDerivRatioFit'});
    
    T_result.Properties.VariableUnits = {'' 'ms' '' '' '' '' 'ms' '' '' '' '' 'ms' '' 'ms' ''};
    
    
    T_result_d = table(y0,peaks_A,peaks_P,...
        sparkParams.amplitude,sparkParams.TTP, sparkParams.TTP_fittedLine, ...
        sparkParams.FDHM,sparkParams.FWHM,sparkParams.sparkMass,...
        sparkParams.tauD,sparkParams.AUC_2DFit,...
        peaksWholeProfFit_A, peaksWholeProfFit_P, peaksWholeProfFit_TTP, peaksWholeProfFit_FDHM,...
        peaksWholeProfFit_AUC ,peaksWholeProfFit_tauD,...
        firstDerMaxVal, firstDerMaxValPos, peaksWholeProfFit_derA, peaksWholeProfFit_derP,...
        'VariableNames',{'baseline' 'Amplitude' 'peakPos'...
        'sparkAmplitude' 'sparkTTP' 'sparkTTP_fittedLine' ...
        'sparkFDHM' 'sparkFWHM' 'sparkMass'...
        'sparkTauD' 'sparkAUC'...
        'AmplitudeFit' 'peakPosFit' 'TTPfit' 'FDHMfit' 'AUCfit' 'tauDfit'...
        'firstDerivMaxVal' 'firstDerivMaxValPos' 'firstDerivMaxValFit' 'firstDerivMaxValPosFit'});
    
    T_result_d.Properties.VariableUnits = {'deltaF/F0' 'deltaF/F0' 'ms'...
        'deltaF/F0' 'ms' 'ms' 'ms' 'um' 'deltaF/F0*um^3' ...
        'ms' 'ms*um*deltaF/F0' ...
        'deltaF/F0' 'ms' 'ms' 'ms' 'ms*deltaF/F0' 'ms' ...
        'deltaF/F0' 'ms' 'deltaF/F0' 'ms'};
    
    result(i,1) = {T_result};
    allPeaksData(i,1) = {T_result_d};
    
    clearvars T_result accepted peaks_A_ratio peaks_A peaks_P peaks_P_diff T_result_d ...
        sparkParams sparkAmp sparkMass sparkMassRatio sparkAmpRatio peaksWholeProfFit_P_diff peaksWholeProfFit_AUC...
        params peaksWholeProfFit_A peaksWholeProfFit_P peaksWholeProfFit_TTP peaksWholeProfFit_FDHM peaksWholeProfFit_tauD ...
        peaksWholeProfFit_A_ratio peaksWholeProfFit_P_ratio peaksWholeProfFit_TTP_ratio peaksWholeProfFit_FDHM_ratio ...
        peaksWholeProfFit_AUC_ratio peaksWholeProfFit_tauD_ratio ...
        peaks_derP_diff peaks_derA_ratio peaksWholeProfFit_derP_diff peaksWholeProfFit_derA_ratio ...
        firstDerMaxVal firstDerMaxValPos peaksWholeProfFit_derA peaksWholeProfFit_derP
    
end

% save data
selectedROIs.AnalysisResult = result;
selectedROIs.allPeaksData = allPeaksData;
profileAnalysis.selectedROIs = selectedROIs;
setappdata(mainFig,'profileAnalysis',profileAnalysis)

% plot result from analysis using fitting only rise part of each
% event
plotProfileAnalysis1(t, selectedROIs, pxSzX)

% plot result from analysis using whole profile fit
if any(strcmp(selectedROIs.Properties.VariableNames, 'wholeProfileFit'))
    plotProfileAnalysis2(t, selectedROIs, pxSzX)
end

% remove empty axes
h_f = findobj('-regexp','Name','repetitive sparks analysis');
h_a = findobj(h_f,'Type','axes');
h_a_d = arrayfun(@(x) isempty(x.Children),h_a);
delete(h_a(h_a_d))

% save data
setappdata(mainFig,'profileAnalysis',profileAnalysis)

% mouse cursor
set(mainFig,'Pointer','arrow')
drawnow

%%%%%%%%%%%%
saveSparkRecAnalysis([],[],mainFig)
%%%%%%%%%%%%

end

