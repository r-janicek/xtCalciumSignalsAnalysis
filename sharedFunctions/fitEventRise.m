function out = fitEventRise(t, prof_t, eventType, options)
%{
parameters of function
t = time vector
prof_t = time profile
eventType = type of event, global or local
options:
    peaks_vals = peaks values
    peaks_locs = peaks positions, in time units 
    ax_prof = handle to axes with time profile
    coefPrevFit = coeficients from previous fit
    fitTol = fit tolerance
    numOfFitIter = number of iterations
    smooth_span = smoothing span prameter
    bs_crit = baseline percentile treshold
    sSpPrev = 
    eSpPrev = 
    evntsMask = mask of events in profile ...
    posOfSelPoints = 
    fitFun = expRise or sigmoid
%}
arguments
    t (:,1) {mustBeNumeric}            
    prof_t (:,1) {mustBeNumeric}             
    eventType (1,1) {mustBeText, ...
        mustBeMember(eventType,{'global','local'})} 
    options.peaks_vals (:,1) {mustBeNumeric} = nan 
    options.peaks_locs (:,1) {mustBeNumeric} = nan
    options.ax_prof = []
    options.coefPrevFit = []
    options.fitTol = 1e-9
    options.numOfFitIter (1,1) {mustBeInteger, mustBePositive} = 1000
    options.smooth_span (1,1) {mustBeInteger, mustBePositive} = 50 % ms
    options.bs_crit (1,1) {mustBeInteger, mustBePositive, ...
        mustBeLessThan(options.bs_crit,100)} = 75 
    options.sSpPrev = []
    options.eSpPrev = []
    options.evntsMask (:,1) logical = true(size(t))
    options.posOfSelPoints = []
    options.fitFun (1,1) {mustBeText, ...
        mustBeMember(options.fitFun,{'expRise','sigmoid'})} = "expRise"
    options.fittedLineColor (:,3) {mustBeNumeric} = [1 0 0]
    options.fittedLineTag (1,1) {mustBeText} = ""
end

% get time step from time vector
pxSzT = mean(diff(t));
if any(isnan(options.peaks_vals)) || any(isnan(options.peaks_locs))
    [options.peaks_vals, options.peaks_locs] = max(prof_t);
    options.peaks_locs = t(options.peaks_locs);
end
% it peaks were selected
if ~isempty(options.posOfSelPoints)
    options.peaks_vals = options.posOfSelPoints.peakFitVal;
    options.peaks_locs = t(options.posOfSelPoints.peakFitPos);
end
% fit options
fitOpts = optimoptions('lsqnonlin', 'TolFun',options.fitTol, ...
    'TolX',options.fitTol, 'MaxIter',options.numOfFitIter,...
    'MaxFunEvals',3*options.numOfFitIter, 'Display','off');
switch options.fitFun
    case 'expRise'
        % parameters [t0 tauR A bs]
        fun_e = @(x,t) ((t>=x(1)).*((1-exp(-(t-x(1))./x(2))).*(x(3)) + x(4)) + ...
            (t<x(1)).*(x(4)));
        fun = @(x,t,ys) ((t>=x(1)).*((1-exp(-(t-x(1))./x(2))).*(x(3)) + x(4)) + ...
            (t<x(1)).*(x(4)))-ys;
        n_coef = 4;
    case 'sigmoid'
        % parameters [bs A m tau] (baseline, amplitude, middle point, time constant)
        fun_e = @(x,t)   x(1) + (x(2)-x(1))./(1+exp(-(t-x(3))./x(4)));
        fun = @(x,t,ys) (x(1) + (x(2)-x(1))./(1+exp(-(t-x(3))./x(4))))-ys;
        n_coef = 4;
end

% preallocate
out.model = options.fitFun;
out.fitFun = fun_e;
out.h_line = zeros(length(options.peaks_vals),1);
out.coef = zeros(length(options.peaks_vals),n_coef);
out.sp_fit = {zeros(length(options.peaks_vals),1)};
out.startOfEvent = zeros(length(options.peaks_vals),1); 
out.endOfEvent = zeros(length(options.peaks_vals),1);
out.detectedEventsMask = false(numel(prof_t),1);
% maximum duration of baseline in points
switch eventType
    case 'global'
        maxDurOfBaseline = ceil(1000/pxSzT);
    case 'local'
        maxDurOfBaseline = ceil(100/pxSzT);
end

% use previous starts and ends of sparks if any, otherwise find them
if ~isempty(options.peaks_vals)
    % analyze and fit all peaks of events
    for i = 1:numel(options.peaks_locs)
        % get position of peak of events in pixels
        [~, peak_loc_px] = min(abs(t-options.peaks_locs(i)));
        if isempty(options.sSpPrev) && isempty(options.eSpPrev)
            % find new start and end of event with baseline
            [pos_s, pos_e] = estimateStartAndEndOfEvent( ...
                prof_t, peak_loc_px, ...
                maxDurOfBaseline=maxDurOfBaseline, ...
                evntsMask=options.evntsMask, ...
                equalBaselineDur=false, ...
                smoothSpan=round(options.smooth_span/pxSzT), ...
                evntAcceptCrit=options.bs_crit);
        else
            % take provided ones
            pos_s = options.sSpPrev(i);
            pos_e = options.eSpPrev(i);
        end
        
        % get rise part of profile of event 
        y_rise = prof_t(pos_s:peak_loc_px);
        t_rise = t(pos_s:peak_loc_px);
        % check if there is enough points to fit the profile with model
        if length(y_rise)<size(out.coef,2)
            y_rise = prof_t(pos_s:peak_loc_px + (size(out.coef,2)-length(y_rise)));
            t_rise = t(pos_s:peak_loc_px + (size(out.coef,2)-length(t_rise)));
        end
        t_rise = t_rise(:);
        y_rise = y_rise(:);

        % get initial parameters to fit profile
        if ~isempty(options.coefPrevFit) && ...
                size(options.coefPrevFit,2) == n_coef
            x0 = options.coefPrevFit(i,:);
        else
            % estimation of initial fit values
            A = options.peaks_vals(i);
            % first estimation of baseline
            try
                bsFirstEstPx = ...
                    find(~options.evntsMask(1:peak_loc_px), 1, 'last') - ...
                    pos_s + 1;
            catch 
                bsFirstEstPx = 3;
            end
            % use first 3 points of profile, if there is a problem or no
            % mask of events
            if isempty(bsFirstEstPx), bsFirstEstPx = 3; end
            if bsFirstEstPx <= 0, bsFirstEstPx = 3; end
            y0 = mean(y_rise(1:bsFirstEstPx));
            try
                % estimate t0
                % find 10 and 90% from max of rise part
                p_90 = numel(y_rise) - ...
                    find( (flipud(y_rise)-(y0 + (max(y_rise)-y0).*0.90))>0, 1, 'last') + 1;
                p_10 = numel(y_rise) - ...
                    find( (flipud(y_rise)-(y0 + (max(y_rise)-y0).*0.10))<0, 1, 'first') + 1;
                %[~,p_75] = min( abs( ys-(y0 + (max(ys)-y0).*0.75) ) );
                %[~,p_25] = min( abs( ys-(y0 + (max(ys)-y0).*0.25) ) );
                if p_90 == p_10
                    p_10 = p_90-1;
                end
                v_90 = y_rise(p_90);
                v_10 = y_rise(p_10);
                % construct line two-point form and get x at y=y0
                t0_est = t_rise(p_10)+( (y0-v_10)*(t_rise(p_90)-t_rise(p_10)) )/(v_90-v_10);
                tR_est = t_rise(p_90)-t_rise(p_10);
                if tR_est<mean(diff(t_rise)), tR_est = 3; end
                if t0_est<t_rise(1), t0_est = t_rise(1); end
            catch
                t0_est = options.peaks_locs(i)-10;
                if t0_est<=0, t0_est = mean(diff(t_rise)); end
                tR_est = 3;
            end
            switch options.fitFun
                case 'expRise'
                    % parameters [t0, tR, A, y0]
                    x0 = [t0_est tR_est A-y0 y0];
                case 'sigmoid'
                    % parameters [bs A m tau]
                    x0 = [y0 A t0_est+(peaks_locs(i)-t0_est)/2 tR_est];
            end
        end
        % parameters bounds
        switch options.fitFun
            case 'expRise'
                % t0, tR, A, y0
                lb = [min(t_rise) 0 min(y_rise) min(y_rise)];
                ub = [max(t_rise) max(t_rise)-min(t_rise) max(y_rise)-min(y_rise) max(y_rise)];
                % lock baseline on selected interval
                if ~isempty(options.posOfSelPoints)
                    selected_bs_val = mean( ...
                        prof_t(options.posOfSelPoints.bsFitS: ...
                               options.posOfSelPoints.bsFitE));
                    x0(4) = selected_bs_val;
                    lb(4) = selected_bs_val;
                    ub(4) = selected_bs_val;
                end
            case 'sigmoid'
                % parameters [bs A m tau]
                lb = [min(y_rise) min(y_rise) min(t_rise) 10^(-9)];
                ub = [max(y_rise) max(y_rise) max(t_rise) inf];
        end
        % fit
        try
            x = lsqnonlin(@(x) fun(x,t_rise,y_rise), x0, lb, ub, fitOpts);
        catch
            try
                sumSqrs = @(x,t,ys) sum(fun(x,t,ys).^2);
                x = fmincon(@(x) sumSqrs(x,t_rise,y_rise), x0, ...
                    [], [], [], [], lb, ub);
            catch
                x = x0;
            end
        end
        % estimate t0 if fit fun was sigmoidal
        switch options.fitFun
            case 'expRise'
            case 'sigmoid'
                % tangent line at middle point of sigmoid
                fun_line = @(x,t) ((x(2)-x(1))/(4*x(4))).*t + ...
                    (fun_e(x,x(3))-(((x(2)-x(1))/(4*x(4))).*x(3)));
                t0 = -(fun_e(x,x(3))-(((x(2)-x(1))/(4*x(4))).*x(3))) / ...
                    ((x(2)-x(1))/(4*x(4)));
        end

        % figure
        % plot(t,ys,'ok')
        % hold on
        % plot(t,fun_e(x0,t),'b')
        % plot(t,fun_e(x,t),'r')
        % %t_ups = linspace(t(1),t(end),1000);
        % %plot(t_ups,fun_e(x,t_ups),'c')
        % 
        % %plot(t(t>t0_est),fun_line(x,t(t>t0_est)),'m')
        
        % save results 
        out.coef(i,:) = x;
        out.sp_fit(i,1) = {[t_rise,fun_e(x,t_rise)]};
        out.startOfEvent(i,1) = pos_s; % in pixels
        out.endOfEvent(i,1) = pos_e; % in pixels
        
        if ~isempty(options.ax_prof)
            out.h_line(i,1) = line(t_rise, fun_e(x,t_rise), ...
                'Parent',options.ax_prof, ...
                'Color',options.fittedLineColor, ...
                'LineStyle','-', ...
                'LineWidth',2, ...
                'Tag',options.fittedLineTag);   
            
            out.detectedEventsMask(pos_s:pos_e,1) = true(pos_e-pos_s+1,1);
        end
    end
    
%     hl_m = line(x_t(detectedEventsMask),prof_t(detectedEventsMask),'Parent',ax_prof,'Color','g',...
%                     'LineStyle','none','Marker','.','MarkerSize',20,'LineWidth',1,'Tag','eventsMask');
%     uistack(hl_m, 'bottom')

end
end


