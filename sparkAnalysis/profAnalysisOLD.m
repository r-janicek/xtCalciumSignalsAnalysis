function profAnalysisOLD(~,~,mainFig)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% get data 
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');
profileAnalysis = getappdata(mainFig,'profileAnalysis');

if ~isfield(profileAnalysis,'selectedROIs')
   return 
end

h = findobj('-regexp','Name','repeated sparks analysis');
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
        
        wholeProfileFit = selectedROIs.wholeProfileFit{i,1};
        % normalized profile
        prof = wholeProfileFit.yN;
        
        locs = cell2mat(wholeProfileFit.eventsPeaks(:,2));
        
        pks = prof(round(locs./pxSzT)+1);
        
        coefPrevFit = selectedROIs.prevFitCoef{i,1};   
        coefPrevFit(:,3) = pks;
        coefPrevFit(:,4) = coefPrevFit(:,4)./mean(wholeProfileFit.y(wholeProfileFit.baselineM));
        
    else
        % normalized profile
        prof = mean(selectedROIs.dataProfile{i,1},1);
        coefPrevFit = selectedROIs.prevFitCoef{i,1};
        locs = cell2mat(selectedROIs.eventsPeaks{i,1}(:,2));    
        pks = prof(round(locs./pxSzT)+1);
    end
            
    dataWholeAreaRaw = selectedROIs.dataWholeAreaRaw{i,1};    
    statEvents = selectedROIs.regionsProp{i,1};   
    posOfROI = selectedROIs.positionOfRoi(i,1);
    smooth_span = selectedROIs.smooth_span(i,1);
    bs_crit = selectedROIs.bs_crit(i,1);
    sSp = selectedROIs.startPosOfSpark{i,1};
    eSp = selectedROIs.endPosOfSpark{i,1};
     
    %fit only rise of spark fun(t0,u,tR,A,y0)    
    [~,~,coef,sp_fit,startOfSpark,endOfSpark] = fitSparkRise(pxSzT,t,prof,pks,locs,[],coefPrevFit,1e-9,1000,smooth_span,bs_crit,sSp,eSp);
    fit_result(i,1) = {coef};
    fit_result(i,2) = {sp_fit};
    
    
    % calculate sparks parameters
    if hObjs.check_sparkParams.Value
        
        % selected method of parameters calculation
        switch hObjs.calcSpParamMethodRBgroup.SelectedObject.String
            
            case '2D gaussian fit'
                calcMethod = '2DGauss';
                
            case 'max crossing profiles'
                calcMethod = 'peakXTProfiles';
        
        end

        % progress bar
        progressbar(sprintf('parameters of events from profile position %d (px)',selectedROIs.positionOfRoi(i)))
        
        %calculate parameters of selected events
        eventsParams = findDetectedEventsParams(dataWholeAreaRaw,statEvents,startOfSpark,endOfSpark,coef,posOfROI,mainFig,calcMethod);
        eventsParameters(i,1) = {eventsParams};
       
    else
        
        eventsParams.amplitude = nan(size(statEvents));
        eventsParams.TTP = nan(size(statEvents));
        eventsParams.FDHM = nan(size(statEvents));
        eventsParams.FWHM = nan(size(statEvents));
        eventsParams.sparkMass = nan(size(statEvents));
        eventsParams.tauD = nan(size(statEvents));
        eventsParams.AUC_2DFit = nan(size(statEvents));
        
        eventsParameters(i,1) = {eventsParams};
        
    end
end

% create result structure
selectedROIs.finalFitResults = fit_result;
selectedROIs.eventsParameters = eventsParameters;


%select type of analysis
if hObjs.check_pairAnalysis.Value    
    pairwise = 'Yes';
else
    pairwise = 'No';    
end

% delays between photolytic pulses
TPP_delays = imgData.TPP_delays;
durTPP = imgData.durOfTPP;


% choose type of analysis of events; pairwise or not 
switch pairwise
    
    case 'Yes'
               
        %analyze peaks ratio and distances
        for i = 1:height(selectedROIs)
            
            if ~isempty(selectedROIs.eventsPeaks{i,1})
                
                % calculate parameters of events from whole profile fit
                if any(strcmp(selectedROIs.Properties.VariableNames,'wholeProfileFit'))
                    paramsFittedEvents = calcParametersOfEventsFromWholeProfileFit(t,selectedROIs.wholeProfileFit{i,1}.profFit);
                    % normalized profiles
                    profN = selectedROIs.wholeProfileFit{i,1}.yN; 
                                        
                else
                    
                    names_params = {'t0 (ms)','tPeak (ms)','Apeak (F/F0)','TTP (ms)','FDHM (ms)','AUC (ms*F/F0)','tauD (ms)','firstDerMaxVal (deltaF/F0)','firstDerMaxValPos (ms)'};
                    paramsFittedEvents = [names_params;num2cell(nan(size(selectedROIs.eventsPeaks{i,1},1),9))]; 
                    % normalized profiles
                    profN = mean(selectedROIs.dataProfile{i,1},1); 
                
                end
                   
                selectedROIs.normProf(i) = {profN};        
                
                keyboard
                % calculate 1. derivative max amplitude and position for
                % all selected events 
                
                
                %detected peak amplitude minus baseline (y0 from fitting of rise of spark)
                y0 = selectedROIs.finalFitResults{i,1}(:,4);
                peaks_P = cell2mat(selectedROIs.eventsPeaks{i,1}(:,2));                
                peaks_A = profN(round(peaks_P./pxSzT)+1) - y0;
                
                % save peaks values
                selectedROIs.eventsPeaks{i,1}(:,1) = num2cell(peaks_A);
                                
                peaksWholeProfFit_A = [paramsFittedEvents{2:end,3}]';
                peaksWholeProfFit_P = [paramsFittedEvents{2:end,2}]'; 
                peaksWholeProfFit_TTP = [paramsFittedEvents{2:end,4}]';
                peaksWholeProfFit_FDHM = [paramsFittedEvents{2:end,5}]';
                peaksWholeProfFit_AUC = [paramsFittedEvents{2:end,6}]';
                peaksWholeProfFit_tauD = [paramsFittedEvents{2:end,7}]';
                
                sparkParams = selectedROIs.eventsParameters{i,1};    
                sparkAmp = sparkParams.amplitude;
                sparkMass = sparkParams.sparkMass;
                
                if ~isempty(peaks_P)
                    
                    for j = 1:numel(peaks_P)/2
               
                        peaks_A_ratio(j,1) = peaks_A(2*j)/peaks_A(2*j-1);
                        peaks_P_diff(j,1) = peaks_P(2*j)-peaks_P(2*j-1);
                        accepted(j,1) = true;
                        sparkMassRatio(j,1) = sparkMass(2*j)/sparkMass(2*j-1);
                        sparkAmpRatio(j,1) = sparkAmp(2*j)/sparkAmp(2*j-1);
                                              
                        peaksWholeProfFit_A_ratio(j,1) = peaksWholeProfFit_A(2*j)/peaksWholeProfFit_A(2*j-1);
                        peaksWholeProfFit_P_diff(j,1) = peaksWholeProfFit_P(2*j)-peaksWholeProfFit_P(2*j-1);
                        peaksWholeProfFit_TTP_ratio(j,1) = peaksWholeProfFit_TTP(2*j)/peaksWholeProfFit_TTP(2*j-1);
                        peaksWholeProfFit_FDHM_ratio(j,1) = peaksWholeProfFit_FDHM(2*j)/peaksWholeProfFit_FDHM(2*j-1);
                        peaksWholeProfFit_AUC_ratio(j,1) = peaksWholeProfFit_AUC(2*j)/peaksWholeProfFit_AUC(2*j-1);
                        peaksWholeProfFit_tauD_ratio(j,1) = peaksWholeProfFit_tauD(2*j)/peaksWholeProfFit_tauD(2*j-1);
                        
                                              
                    end
                    
                else
                    peaks_A_ratio = [];
                    peaks_P_diff = [];
                    accepted = [];
                    sparkMassRatio = [];
                    sparkAmpRatio = [];
                    
                    peaksWholeProfFit_A_ratio = [];
                    peaksWholeProfFit_P_diff = [];
                    peaksWholeProfFit_TTP_ratio = [];
                    peaksWholeProfFit_FDHM_ratio = [];
                    peaksWholeProfFit_AUC_ratio = [];
                    peaksWholeProfFit_tauD_ratio = [];
                    
                end
                
            else
                peaks_A_ratio = [];
                peaks_P_diff = [];
                accepted = [];
                sparkMassRatio = [];
                sparkAmpRatio = [];
                
                peaksWholeProfFit_A_ratio = [];
                peaksWholeProfFit_P_diff = [];
                peaksWholeProfFit_TTP_ratio = [];
                peaksWholeProfFit_FDHM_ratio = [];
                peaksWholeProfFit_AUC_ratio = [];
                peaksWholeProfFit_tauD_ratio = [];

            end
      
            T_result = table(peaks_A_ratio,peaks_P_diff,accepted,sparkAmpRatio,sparkMassRatio,...
                             peaksWholeProfFit_A_ratio,peaksWholeProfFit_P_diff,peaksWholeProfFit_TTP_ratio,...
                             peaksWholeProfFit_FDHM_ratio,peaksWholeProfFit_AUC_ratio,peaksWholeProfFit_tauD_ratio,...
                'VariableNames',{'AmplitudeRatio' 'peakPosDiff' 'acceptedPair' 'sparkAmpRatio' 'sparkMassRatio' ... 
                                 'AmplRatioFit' 'peakPosDiffFit' 'TTPratioFit' 'FDHMratioFit' 'AUCratio' 'tauDratioFit'});
            
            T_result.Properties.VariableUnits = {'' 'ms' '' '' '' '' 'ms' '' '' '' ''};
     
            
            T_result_d = table(y0,peaks_A,peaks_P,...
                              sparkParams.amplitude,sparkParams.TTP,sparkParams.FDHM,sparkParams.FWHM,sparkParams.sparkMass,...
                              sparkParams.tauD,sparkParams.AUC_2DFit,...
                              peaksWholeProfFit_A, peaksWholeProfFit_P, peaksWholeProfFit_TTP, peaksWholeProfFit_FDHM,...
                              peaksWholeProfFit_AUC,peaksWholeProfFit_tauD,...
                              'VariableNames',{'baseline' 'Amplitude' 'peakPos'...
                                               'sparkAmplitude' 'sparkTTP' 'sparkFDHM' 'sparkFWHM' 'sparkMass' ...
                                               'sparkTauD' 'sparkAUC'...
                                               'AmplitudeFit' 'peakPosFit' 'TTPfit' 'FDHMfit' 'AUCfit' 'tauDfit'});
                        
            T_result_d.Properties.VariableUnits = {'deltaF/F0' 'deltaF/F0' 'ms'...
                                                   'deltaF/F0' 'ms' 'ms' 'um' 'deltaF/F0*um^3' ...
                                                   'ms' 'ms*um*deltaF/F0' ...
                                                   'deltaF/F0' 'ms' 'ms' 'ms' 'ms*deltaF/F0' 'ms'};
            
            result(i,1) = {T_result};
            allPeaksData(i,1) = {T_result_d};
            
            clearvars T_result accepted peaks_A_ratio peaks_A peaks_P peaks_P_diff T_result_d ...
                      sparkParams sparkAmp sparkMass sparkMassRatio sparkAmpRatio ...
                      params peaksWholeProfFit_A peaksWholeProfFit_P peaksWholeProfFit_TTP ...
                      peaksWholeProfFit_FDHM peaksWholeProfFit_AUC peaksWholeProfFit_tauD ...
                      peaksWholeProfFit_A_ratio peaksWholeProfFit_P_ratio peaksWholeProfFit_TTP_ratio ...
                      peaksWholeProfFit_FDHM_ratio peaksWholeProfFit_tauD_ratio peaksWholeProfFit_P_diff peaksWholeProfFit_AUC_ratio
                  
        end
          
        % save data
        selectedROIs.AnalysisResult = result;
        selectedROIs.allPeaksData = allPeaksData;
        
        profileAnalysis.selectedROIs = selectedROIs;
        
        setappdata(mainFig,'profileAnalysis',profileAnalysis)
               
        % plot result from analysis using fitting only rise part of each
        % event
        plotProfileAnalysis1(t,selectedROIs,TPP_delays,durTPP,pxSzX,pairwise,getappdata(mainFig,'analysisType'))
        
        % plot result from analysis using whole profile fit
        if any(strcmp(selectedROIs.Properties.VariableNames,'wholeProfileFit'))
            plotProfileAnalysis2(t,selectedROIs,TPP_delays,durTPP,pxSzX,pairwise,getappdata(mainFig,'analysisType'))
        end
        
    case 'No'
                
        %analyze peaks ratio and distances
        for i = 1:height(selectedROIs)
            
            if ~isempty(selectedROIs.eventsPeaks{i,1})
     
                 % calculate parameters of events from whole profile fit
                if any(strcmp(selectedROIs.Properties.VariableNames,'wholeProfileFit'))
                    paramsFittedEvents = calcParametersOfEventsFromWholeProfileFit(t,selectedROIs.wholeProfileFit{i,1}.profFit);
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
                
                keyboard
                % calculate 1. derivative max amplitude and position for
                % all selected events 
                figure
                plot(profN,'k')
                hold on 
                plot(gradient(profN),'r')
                
                
                %detected peak amplitude minus baseline (y0 from fitting of rise of spark)
                y0 = selectedROIs.finalFitResults{i,1}(:,4);
                peaks_P = cell2mat(selectedROIs.eventsPeaks{i,1}(:,2));                
                peaks_A = profN(round(peaks_P./pxSzT)+1)-y0;
                
                % save peaks values
                selectedROIs.eventsPeaks{i,1}(:,1) = num2cell(peaks_A);
                
                % peaks positions difference             
                peaks_P_diff = diff(peaks_P);
        
                sparkParams = [selectedROIs.eventsParameters{i,1}];    
                sparkAmp = sparkParams.amplitude;
                sparkMass = sparkParams.sparkMass;
           
                peaksWholeProfFit_A = [paramsFittedEvents{2:end,3}]';
                peaksWholeProfFit_P = [paramsFittedEvents{2:end,2}]'; 
                peaksWholeProfFit_TTP = [paramsFittedEvents{2:end,4}]';
                peaksWholeProfFit_FDHM = [paramsFittedEvents{2:end,5}]';
                peaksWholeProfFit_AUC = [paramsFittedEvents{2:end,6}]';
                peaksWholeProfFit_tauD = [paramsFittedEvents{2:end,7}]';
                              
                if ~isempty(peaks_P_diff)
                    
                    for j = 1:length(peaks_P_diff)
                        
                        peaks_A_ratio(j,1) = peaks_A(j+1)/peaks_A(j);
                                              
                        sparkMassRatio(j,1) = sparkMass(j+1)/sparkMass(j);
                        sparkAmpRatio(j,1) = sparkAmp(j+1)/sparkAmp(j);
                                               
                        peaksWholeProfFit_A_ratio(j,1) = peaksWholeProfFit_A(j+1)/peaksWholeProfFit_A(j);
                        peaksWholeProfFit_P_diff(j,1) = peaksWholeProfFit_P(j+1)-peaksWholeProfFit_P(j);
                        peaksWholeProfFit_TTP_ratio(j,1) = peaksWholeProfFit_TTP(j+1)/peaksWholeProfFit_TTP(j);
                        peaksWholeProfFit_FDHM_ratio(j,1) = peaksWholeProfFit_FDHM(j+1)/peaksWholeProfFit_FDHM(j);
                        peaksWholeProfFit_AUC_ratio(j,1) = peaksWholeProfFit_AUC(j+1)/peaksWholeProfFit_AUC(j);
                        peaksWholeProfFit_tauD_ratio(j,1) = peaksWholeProfFit_tauD(j+1)/peaksWholeProfFit_tauD(j);
                        
                                                
                        if j==1
                     
                            if peaks_P(1) < 200
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
                    accepted = [];                    
                    sparkMassRatio = [];
                    sparkAmpRatio = [];
                    
                    peaksWholeProfFit_A_ratio = [];
                    peaksWholeProfFit_P_diff = [];
                    peaksWholeProfFit_TTP_ratio = [];
                    peaksWholeProfFit_FDHM_ratio = [];
                    peaksWholeProfFit_AUC_ratio = [];
                    peaksWholeProfFit_tauD_ratio = [];
                    
                end
                
            else
                peaks_A_ratio = [];
                peaks_P_diff = [];
                accepted = [];
                sparkMassRatio = [];
                sparkAmpRatio = [];
                
                peaksWholeProfFit_A_ratio = [];
                peaksWholeProfFit_P_diff = [];
                peaksWholeProfFit_TTP_ratio = [];
                peaksWholeProfFit_FDHM_ratio = [];
                peaksWholeProfFit_AUC_ratio = [];
                peaksWholeProfFit_tauD_ratio = [];
                
            end
                                                          
            T_result = table(peaks_A_ratio,peaks_P_diff,accepted,sparkAmpRatio,sparkMassRatio,...
                             peaksWholeProfFit_A_ratio,peaksWholeProfFit_P_diff,peaksWholeProfFit_TTP_ratio,...
                             peaksWholeProfFit_FDHM_ratio,peaksWholeProfFit_AUC_ratio,peaksWholeProfFit_tauD_ratio,...
                'VariableNames',{'AmplitudeRatio' 'peakPosDiff' 'acceptedPair' 'sparkAmpRatio' 'sparkMassRatio' ... 
                                 'AmplRatioFit' 'peakPosDiffFit' 'TTPratioFit' 'FDHMratioFit' 'AUCratio' 'tauDratioFit'});
            
            T_result.Properties.VariableUnits = {'' 'ms' '' '' '' '' 'ms' '' '' '' ''};
            
         
            T_result_d = table(y0,peaks_A,peaks_P,...
                              sparkParams.amplitude,sparkParams.FDHM,sparkParams.TTP,sparkParams.FWHM,sparkParams.sparkMass,...
                              sparkParams.tauD,sparkParams.AUC_2DFit,...
                              peaksWholeProfFit_A, peaksWholeProfFit_P, peaksWholeProfFit_TTP, peaksWholeProfFit_FDHM,...
                              peaksWholeProfFit_AUC,peaksWholeProfFit_tauD,...
                              'VariableNames',{'baseline' 'Amplitude' 'peakPos'...
                                               'sparkAmplitude' 'sparkFDHM' 'sparkTTP' 'sparkFWHM' 'sparkMass'...
                                               'sparkTauD' 'sparkAUC'...
                                               'AmplitudeFit' 'peakPosFit' 'TTPfit' 'FDHMfit' 'AUCfit' 'tauDfit'});
                         
            T_result_d.Properties.VariableUnits = {'deltaF/F0' 'deltaF/F0' 'ms'...
                                                   'deltaF/F0' 'ms' 'ms' 'um' 'deltaF/F0*um^3' ...
                                                   'ms' 'ms*um*deltaF/F0' ...
                                                   'deltaF/F0' 'ms' 'ms' 'ms' 'ms*deltaF/F0' 'ms'};
            
            result(i,1) = {T_result};
            allPeaksData(i,1) = {T_result_d};
            
            clearvars T_result accepted peaks_A_ratio peaks_A peaks_P peaks_P_diff T_result_d ...
                      sparkParams sparkAmp sparkMass sparkMassRatio sparkAmpRatio peaksWholeProfFit_P_diff peaksWholeProfFit_AUC...
                      params peaksWholeProfFit_A peaksWholeProfFit_P peaksWholeProfFit_TTP peaksWholeProfFit_FDHM peaksWholeProfFit_tauD ...
                      peaksWholeProfFit_A_ratio peaksWholeProfFit_P_ratio peaksWholeProfFit_TTP_ratio peaksWholeProfFit_FDHM_ratio ...
                      peaksWholeProfFit_AUC_ratio peaksWholeProfFit_tauD_ratio
                                                                                                                   
        end
               
        % save data
        selectedROIs.AnalysisResult = result;
        selectedROIs.allPeaksData = allPeaksData;
        
        profileAnalysis.selectedROIs = selectedROIs;        
        setappdata(mainFig,'profileAnalysis',profileAnalysis)
          
        % plot result from analysis using fitting only rise part of each
        % event
        plotProfileAnalysis1(t,selectedROIs,TPP_delays,durTPP,pxSzX,pairwise,getappdata(mainFig,'analysisType'))
        
        % plot result from analysis using whole profile fit
        if any(strcmp(selectedROIs.Properties.VariableNames,'wholeProfileFit'))
            plotProfileAnalysis2(t,selectedROIs,TPP_delays,durTPP,pxSzX,pairwise,getappdata(mainFig,'analysisType'))
        end
end

h_f = findobj('-regexp','Name','repeated sparks analysis');
h_a = findobj(h_f,'Type','axes');
h_a_d = arrayfun(@(x) isempty(x.Children),h_a);
delete(h_a(h_a_d))

% save data
profileAnalysis.pairwise = pairwise;
profileAnalysis.TPP_delays = TPP_delays;

setappdata(mainFig,'profileAnalysis',profileAnalysis)

% mouse cursor
set(mainFig,'Pointer','arrow')
drawnow


%%%%%%%%%%%%
saveSparkRecAnalysis([],[],mainFig)
%%%%%%%%%%%%


end

