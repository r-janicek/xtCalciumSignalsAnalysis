function params = calcParametersOfEventsFromWholeProfileFit(x_t,fitData)
% profile is normalized
% calculate parameters:
% [t0 (start of event), tPeak (time pos of peak), 
% Apeak (amplitude of peak), TTP (time to peak), 
% FDHM (full duration in half maximum)] of events from whole profile fit 
if isfield(fitData,'t_ups')
    % upscaled time axes
    allEventsFits = cell2mat(fitData.t_ups.individualEventsFits(2:end,:));
    t = fitData.t_ups.t_ups;
else
    allEventsFits = cell2mat(fitData.individualEventsFits(2:end,:));
    t = x_t;  
end

% exp decay function 
funExpFit = @(p,x,y) p(1)*exp(-(x-p(2))/p(3)) - y;
funExpEval = @(p,x) p(1)*exp(-(x-p(2))/p(3));

eventsParams = zeros(size(allEventsFits,2),9);
for i = 1:size(allEventsFits,2)
    
    % single event
    eventFit = allEventsFits(:,i);
    % calculate time parameters of single event
    m_coeff_t0 = strcmp(fitData.coefficientsOfFittedEvents(1,:), 't0');
    if any(m_coeff_t0)
        t0 = fitData.coefficientsOfFittedEvents{i+1,m_coeff_t0};
    else
        t0 = [];
    end
    eventParams = getParametersOfEventProfile( ... 
        t, eventFit, [], [], 0, 0, t0, 0, [], [], 1);
    
    % calculate area under curve
    AUC = trapz(t, eventFit);
    
    % calculate 1. derivative max amplitude and position
    dt = mean(diff(t));
    [firstDerVal, firstDerPos] = max(gradient(eventFit)./dt); 
       
    % get tau of decay
    try
        ind = strcmp(fitData.coefficientsOfFittedEvents(1,:),'tauD');
        try 
            tauD = fitData.coefficientsOfFittedEvents{i+1,ind};
        catch 
            tauD = [];
        end
        
        ind_r = strcmp(fitData.coefficientsOfFittedEvents(1,:),'tauR');
        try 
            tauR = fitData.coefficientsOfFittedEvents{i+1,ind_r};
        catch 
            tauR = [];
        end
        
        if ~isempty(tauR) && ~isempty(tauD) && tauR > tauD
            tauD = [];
        end
        
        if isempty(tauD)
            % fit with single exponential             
           t_exp = fitData.t;
           y_exp = cell2mat(fitData.individualEventsFits(2:end,i));
           [vM,pM] = max(y_exp);
       
           coefExp = lsqnonlin(@(p)funExpFit(p,t_exp(pM:end), ...
               y_exp(pM:end)), ...
               [vM t_exp(pM) 20], ...
               [0 0 0], ...
               [inf inf inf]);
           
           tauD = coefExp(3);          
        end               
    catch
        tauD = nan;                       
    end
  
    % final matrix of parameters
    try
        eventsParams(i,:) = [eventParams.t0, eventParams.t_max, ...
                             eventParams.v_max, eventParams.TTP, ...
                             eventParams.FDHM, AUC, ...
                             tauD, firstDerVal, t(firstDerPos)];
    catch
        eventsParams(i,:) = [nan, nan, nan, nan, nan, nan, nan, nan, nan];
    end  
end

% save data 
names_params = {'t0 (ms)', 'tPeak (ms)', 'Apeak (deltaF/F0)', ...
                'TTP (ms)', 'FDHM (ms)', 'AUC (ms*deltaF/F0)', ...
                'tauD (ms)', 'firstDerMaxVal (deltaF/F0)', ...
                'firstDerMaxValPos (ms)'};
params = [names_params;num2cell(eventsParams)];

end

