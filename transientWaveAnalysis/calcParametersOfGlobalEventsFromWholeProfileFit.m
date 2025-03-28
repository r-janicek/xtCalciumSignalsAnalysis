function params = calcParametersOfGlobalEventsFromWholeProfileFit...
    (t, splFit, t_ups, splFit_ups, coefFitEventRise, ...
    startOfEvent, endOfEvent, Npeaks, pks, locs, expFitS, expFitE)

%{ 
t = real time axes
t_ups = upscaled time axes from spline fit
splFit = upscaled spline fit
coefFitEventRise = coef from previous fitting [t0, tauR, A, bs]
startOfEvent =  in pixels
endOfEvent = in pixels
expFitS = in points
expFitE = in points

calculate parameters [t0 (start of event), tPeak (time pos of peak), 
Apeak (amplitude of peak), TTP (time to peak), 
FDHM (full duration in half maximum)] of events from whole profile fit 

%}

% params from event rise fit 
bs = coefFitEventRise(:,4);
t0 = coefFitEventRise(:,1);
tauR = coefFitEventRise(:,2);

t = t(:);
if isempty(t_ups)
    t_ups = t;
end   
t_ups = t_ups(:);
splFit = splFit(:);
splFit_ups = splFit_ups(:);

t_indx = (1:1:length(t_ups));
t_indx = t_indx(:);

% exp decay function 
funExpFit = @(p,x,y) (p(1)*exp(-(x-p(2))/p(3))+p(4)) - y;
funExpEval = @(p,x) p(1)*exp(-(x-p(2))/p(3))+p(4);

eventsParams = zeros(Npeaks, 17);
for i = 1:Npeaks
    try 
        % start of event in points (in t_ups)
        [~, eS] = min( abs( t_ups - t(startOfEvent(i)) ) ); 
        % end of event in points (in t_ups)
        %[~, eE] = min( abs( t_ups - t(endOfEvent(i)) ) );
        [~, eE] = min( abs( t_ups - (t0(i+1) - 50) ) ); % 50 ms before start of event

        % in real time axis, in points
        eSr = startOfEvent(i);
        % eEr = endOfEvent(i);
        [~,eEr] = min( abs( t - (t0(i+1) - 50) ) ); % 50 ms before start of event
    catch
        % start of event in points (in t_ups)
        [~,eS] = min( abs( t_ups - t(startOfEvent(i)) ) ); 
        % end of event in points (in t_ups)
        eE = numel(t_ups);
        
        % in real time axis, in points
        eSr = startOfEvent(i);
        eEr = numel(t); 
    end
    
    if eS<2, eS = 2; end
    if eSr<2, eSr = 2; end
    
    % single event with same length as t_ups or t
    eventFit_ups = [ones(eS-1,1).*bs(i); 
                    splFit_ups(eS:eE); 
                    ones(numel(t_ups)-eE,1).*bs(i) ];
    eventFit = [ones(eSr-1,1).*bs(i); 
                splFit(eSr:eEr); 
                ones(numel(t)-eEr,1).*bs(i) ];
    % position and value of peak
    [pks_ups,locs_ups] = findpeaks(eventFit_ups);
    % take peak closest to manually chosen one
    [~,idx_locs_in_t_ups] = min(abs(t_ups-locs(i)));

    [~,idx_pM] = min(abs(locs_ups-idx_locs_in_t_ups));
    %[vM,pM] = max(eventFit_ups);
    pM = locs_ups(idx_pM);
    vM = pks_ups(idx_pM);
    peakData.val = vM;
    peakData.pos = pM;
    % calculate time params of single event
    eP = getParametersOfEventProfile( ... 
        t_ups, eventFit_ups, [], [], ...
        bs(i), [], t0(i), 0, peakData, [], 1);

    % calculate area under curve, subtract baseline integral (basically rectangle) 
    AUC = trapz(t_ups, eventFit_ups-bs(i));
    
    % calculate 1. derivative max amplitude and position, in real time axis   
    [firstDerVal, firstDerPos] = max(gradient(eventFit));              
    % get tau of decay
    try
        % fit with single exponential
        if ~isempty(expFitS)
            [vMr, pMr] = max(eventFit);
            t_exp = t(expFitS:expFitE);
            y_exp = eventFit(expFitS:expFitE);
            posExp = 1;
        else
            [vMr,pMr] = max(eventFit);
            t_exp = t(pMr:eEr);
            y_exp = eventFit(pMr:eEr);
            [~,posExp] = min(abs(y_exp - (bs(i)+(vMr-bs(i))*0.85)));
        end
        [~,val_33] = min(abs(y_exp - (bs(i)+(vMr-bs(i))*0.33)));
        tauD_est = t_exp(val_33)-t_exp(1);
        
        opt = optimoptions('lsqnonlin', 'MaxIterations',1000, ...
            'MaxFunctionEvaluations',3000,'FunctionTolerance',1e-9, ...
            'TolFun',1e-9, 'TolX',1e-9);
        coefExp = lsqnonlin(@(p)funExpFit(p,t_exp(posExp:end), ...
            y_exp(posExp:end)),...
            [vMr, t_exp(posExp), tauD_est, bs(i)], ...
            [0 0 0 -inf], ...
            [inf inf inf inf], ...
            opt);
        tauD = coefExp(3);
        
        % figure
        % plot(t_exp,y_exp,'ok')
        % hold on
        % plot(t_exp,funExpEval(coefExp,t_exp),'r')
        expDecayFit(i,1) = {[t_exp,funExpEval(coefExp,t_exp)]}; 
    catch
        tauD = nan;
        expDecayFit(i,1) = {nan};
    end

    % final matrix of parameters
    eP_fiedN = fieldnames(eP);
    try
        eventsParams(i,:) = [
            eP.t0, eP.(eP_fiedN{contains(eP_fiedN, 't0_lineFit_')}), ...
            eP.bs_t, eP.t_max, ... 
            eP.half_max_t_1, eP.half_max_t_2, eP.half_max_t, ...
            eP.v_max, eP.amplitude, ... 
            eP.TTP, eP.(eP_fiedN{contains(eP_fiedN, 'TTP_lineFit_')}), ...
            eP.FDHM, AUC, tauR(i), tauD, firstDerVal, t(firstDerPos)];
    catch
        eventsParams(i,:) = [nan, nan, ...
                             nan, nan, ...
                             nan, nan, nan, ...
                             nan, nan, ...
                             nan, nan, ...
                             nan, nan, nan, nan, nan, nan];
    end
    
    clearvars y_t1 y_t2 eP
end

% save data 
names_params = { ...
    't0 (ms)', [eP_fiedN{contains(eP_fiedN, 't0_lineFit_')},' (ms)'], ...
    'F0 (deltaF/F0)', 'tPeak (ms)', ...
    'halfMax1 (ms)', 'halfMax2 (ms)', 'halfMax (deltaF/F0)',...
    'Apeak (deltaF/F0)', 'amplitude (deltaF/F0)', ...
    'TTP (ms)', [eP_fiedN{contains(eP_fiedN, 'TTP_lineFit_')},' (ms)'] ...
    'FDHM (ms)', 'AUC (ms*deltaF/F0)', 'tauR (ms)', 'tauD (ms)', ...
    'firstDerMaxVal (deltaF/F0)', 'firstDerMaxValPos (ms)', ...
    'expDecayFit'};
params = [names_params; [num2cell(eventsParams), expDecayFit]];

end

