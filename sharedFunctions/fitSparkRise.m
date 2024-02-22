function [h_line, detectedEventsMask, coef, sp_fit, ...
    startOfSpark, endOfSpark] = fitSparkRise(pxSz_t, x_t, prof_t, ...
    peaks_vals, peaks_locs, ax_prof, coefPrevFit, tol, iter,...
    smooth_span, bs_crit, sSpPrev, eSpPrev, prof_t_evnts_m)

% locs in time units
options = optimoptions('lsqnonlin', 'TolFun',tol, ...
    'TolX',tol, 'MaxIter',iter,...
    'MaxFunEvals',3*iter, 'Display','off');
fitFun = 'expRise';
switch fitFun
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

% allocate
h_line = zeros(length(peaks_vals),1);
coef = zeros(length(peaks_vals),n_coef); 
sp_fit = {zeros(length(peaks_vals),1)};
startOfSpark = zeros(length(peaks_vals),1); 
endOfSpark = zeros(length(peaks_vals),1);
detectedEventsMask = false(numel(prof_t),1);

maxDurOfBaseline = ceil(100/pxSz_t); % maximum duration of baseline in points

% use previous starts and ends of sparks if any, otherwise find them
if ~isempty(peaks_vals)
    % analyze and fit all peaks of events
    for i = 1:numel(peaks_locs)
        % get position of peak of events in pixels
        [~,peak_loc_px] = min(abs(x_t-peaks_locs(i)));
        if isempty(sSpPrev) && isempty(eSpPrev)
            % find new start and end of event with baseline
            [pos_s, pos_e] = estimateStartAndEndOfEvent( ...
                prof_t, peak_loc_px, ...
                maxDurOfBaseline=maxDurOfBaseline, ...
                evntsMask=prof_t_evnts_m, ...
                equalBaselineDur=false, ...
                smoothSpan=round(smooth_span/pxSz_t), ...
                evntAcceptCrit=bs_crit);
        else
            % take provided ones
            pos_s = sSpPrev(i);
            pos_e = eSpPrev(i);
        end
        
        % get rise part of profile of event 
        ys = prof_t(pos_s:peak_loc_px);
        t = x_t(pos_s:peak_loc_px);
        % check if there is enough points to fit the profile with model
        if length(ys)<size(coef,2)
            ys = prof_t(pos_s:peak_loc_px + (size(coef,2)-length(ys)));
            t = x_t(pos_s:peak_loc_px + (size(coef,2)-length(t)));
        end
        t = t(:);
        ys = ys(:);

        % get initial parameters to fit profile
        if ~isempty(coefPrevFit) && size(coefPrevFit,2) == n_coef
            x0 = coefPrevFit(i,:);
        else
            % estimation of initial fit values
            A = peaks_vals(i); 
            % t0FirstEst = locs(i)-10;
            % [~,t0FirstEstPx] = min( abs( t-t0FirstEst ) );
            t0FirstEstPx = find(~prof_t_evnts_m(1:peak_loc_px), 1, 'last') - ...
                pos_s + 1;
            if t0FirstEstPx <= 0
                t0FirstEstPx = 3; % first 3 points
            end
            y0 = mean(ys(1:t0FirstEstPx));
            try
                % estimate t0
                % find 25 and 75% from max of rise part
                p_75 = numel(ys) - ...
                    find( (flipud(ys)-(y0 + (max(ys)-y0).*0.75))>0, 1, 'last') + 1;
                p_25 = numel(ys) - ...
                    find( (flipud(ys)-(y0 + (max(ys)-y0).*0.25))<0, 1, 'first') + 1;
                %[~,p_75] = min( abs( ys-(y0 + (max(ys)-y0).*0.75) ) );
                %[~,p_25] = min( abs( ys-(y0 + (max(ys)-y0).*0.25) ) );
                if p_75 == p_25
                    p_25 = p_75-1;
                end
                v_75 = ys(p_75);
                v_25 = ys(p_25);
                % construct line two-point form and get x at y=y0
                t0_est = t(p_25)+( (y0-v_25)*(t(p_75)-t(p_25)) )/(v_75-v_25);
                tR_est = t(p_75)-t(p_25);
                if tR_est<mean(diff(t)), tR_est = 3; end
                if t0_est<t(1), t0_est = t(1); end
            catch
                t0_est = peaks_locs(i)-10;
                tR_est = 3;
            end
            switch fitFun
                case 'expRise'
                    % parameters [t0, tR, A, y0]
                    x0 = [t0_est tR_est A-y0 y0];
                case 'sigmoid'
                    % parameters [bs A m tau]
                    x0 = [y0 A t0_est+(peaks_locs(i)-t0_est)/2 tR_est];
            end
        end
        % parameters bounds
        switch fitFun
            case 'expRise'
                % t0, tR, A, y0
                lb = [min(t) 0 min(ys) min(ys)];
                ub = [max(t) max(t)-min(t) max(ys)-min(ys) max(ys)];
            case 'sigmoid'
                % parameters [bs A m tau]
                lb = [min(ys) min(ys) min(t) 10^(-9)];
                ub = [max(ys) max(ys) max(t) inf];
        end
        % fit
        try
            x = lsqnonlin(@(x) fun(x,t,ys), x0, lb, ub, options);
        catch
            try
                sumSqrs = @(x,t,ys) sum(fun(x,t,ys).^2);
                x = fmincon(@(x) sumSqrs(x,t,ys), x0, ...
                    [], [], [], [], lb, ub);
            catch
                x = x0;
            end
        end
        % estimate t0 if fit fun was sigmoidal
        switch fitFun
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

        coef(i,:) = x;
        sp_fit(i,1) = {[t,fun_e(x,t)]};
        startOfSpark(i,1) = pos_s; % in pixels
        endOfSpark(i,1) = pos_e; % in pixels
        
        if ~isempty(ax_prof)
            h_line(i,1) = line(t,fun_e(x,t),'Parent',ax_prof,'Color','r','LineStyle','-','LineWidth',2);      
            
            %h_lineDetectedEvent(i,1) = line(x_t(pos_s:pos_e),prof_t(pos_s:pos_e),'Parent',ax_prof,'Color','g','LineStyle','--','LineWidth',1);
            detectedEventsMask(pos_s:pos_e,1) = true(pos_e-pos_s+1,1);                          
        end
        
    end
    
%     hl_m = line(x_t(detectedEventsMask),prof_t(detectedEventsMask),'Parent',ax_prof,'Color','g',...
%                     'LineStyle','none','Marker','.','MarkerSize',20,'LineWidth',1,'Tag','eventsMask');
%     uistack(hl_m, 'bottom')
    
    
end

end












