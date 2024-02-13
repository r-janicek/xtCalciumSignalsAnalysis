function [h_line, detectedEventsMask, coef, sp_fit, ...
    startOfSpark, endOfSpark] = fitSparkRise(pxSz_t, x_t, prof_t, ...
    pks, locs, ax_prof, coefPrevFit, tol, iter,...
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
h_line = zeros(length(pks),1);
coef = zeros(length(pks),n_coef); 
sp_fit = {zeros(length(pks),1)};
startOfSpark = zeros(length(pks),1); 
endOfSpark = zeros(length(pks),1);
detectedEventsMask = false(numel(prof_t),1);

maxDurOfBaseline = round(100/pxSz_t); % maximum duration of baseline in points

% use previous starts and ends of sparks if any, otherwise find them
if ~isempty(pks)
    % analyze and fit all peaks of events
    for i = 1:numel(locs)
        % get position of peak of events in pixels
        [~,locs_px] = min(abs(x_t-locs(i)));
        if isempty(sSpPrev) && isempty(eSpPrev)
            % mask of fitted event
            m_event = false(size(prof_t));
            m_event(find(~prof_t_evnts_m(1:locs_px),1,'last'): ...
                locs_px+find(~prof_t_evnts_m(locs_px:end),1,'first')-1) = true;
            % take event peak position and its surrounding
            % (maxDurOfBaseline, 2*maxDurOfBaseline)
            m_eventWithBsl = false(size(prof_t));
            m_eventWithBsl( ...
                max(1, locs_px - maxDurOfBaseline) : ...
                min(numel(prof_t), locs_px + 2*maxDurOfBaseline)) = true;
            % smooth profile for further analysis,
            % loess with defined duration in ms
            %prof_s = smooth(prof_t,3);
            n_pts = round(smooth_span/pxSz_t);
            % remove baseline, only events bigger than specified percentile stay
            prof_s = nan(size(prof_t));
            prof_s(m_eventWithBsl) = smooth(prof_t(m_eventWithBsl), ...
                n_pts/numel(prof_t(m_eventWithBsl)), 'loess');
            bs_crit = round(bs_crit);
            % calculate from event
            percentl = prctile(prof_s(m_eventWithBsl & prof_s>0), ...
                [25 50 bs_crit]);
            if max(percentl) > max(prof_s(prof_t_evnts_m))
                percentl = prctile(prof_s(prof_t_evnts_m), ...
                    [25 50 bs_crit]);
            end
            %iqr = percentl(3)- percentl(1);
            bsl = percentl(3);
            prof_s(isnan(prof_s)) = bsl;
            % treshold profile
            prof_s(prof_s < bsl) = bsl;
            % get all posible peaks in smoothed profile of events
            [valPeaks_s, locPeaks_s] = findpeaks(prof_s(m_event));
            locPeaks_s = locPeaks_s + find(m_event, 1, 'first') - 1;
            % get the one closest to peak of event
            [~, idx_p] = min(abs(locPeaks_s-locs_px));
            % change values of currently fitted event to peak value 
            % (to be sure I have correct estimate of baseline. 
            % It might happen that detected spark has multiple peaks, 
            % this will remove them and still keep nice trace to calculate gradient)
            if numel(locPeaks_s) > 1
                prof_s(m_event) = valPeaks_s(idx_p);
            end
            % calculate start and end of event using derivation
            % calculated on thresholded event profile
            % find start of event
            prof_s_beforePeak = prof_s(1:locPeaks_s(idx_p));
            % flip profile to calulate gradient from left to right
            prof_s_beforePeak = flipud(prof_s_beforePeak(:));
            pos_s = max( [...
                locPeaks_s(idx_p)-find(gradient(prof_s_beforePeak)>0, 1, 'first')+1, ...
                locPeaks_s(idx_p)-maxDurOfBaseline ] );
            if isempty(pos_s) || isnan(pos_s) || pos_s < 1
                pos_s = find(m_eventWithBsl, 1, 'first');
                if isempty(pos_s)
                    pos_s = 1;
                end
            end
            % find end of event
            prof_s_afterPeak = prof_s(locPeaks_s(idx_p):end);
            pos_e = min( [...
                locPeaks_s(idx_p)+find(gradient(prof_s_afterPeak)>0, 1, 'first')-1, ...
                locPeaks_s(idx_p)+2*maxDurOfBaseline] );
            if isempty(pos_e) || isnan(pos_e) || pos_e>numel(prof_t)
                pos_e = find(m_eventWithBsl, 1, 'last');
            end
            if pos_e<=pos_s, pos_e = pos_s+1; end

        else
            % take old
            pos_s = sSpPrev(i);
            pos_e = eSpPrev(i);
        end
        
        % get rise part of profile of event 
        ys = prof_t(pos_s:locs_px);
        t = x_t(pos_s:locs_px);
        % check if there is enough points to fit the profile with model
        if length(ys)<size(coef,2)
            ys = prof_t(pos_s:locs_px + (size(coef,2)-length(ys)));
            t = x_t(pos_s:locs_px + (size(coef,2)-length(t)));
        end
        t = t(:);
        ys = ys(:);

        % get initial parameters to fit profile
        if ~isempty(coefPrevFit) && size(coefPrevFit,2) == n_coef
            x0 = coefPrevFit(i,:);
        else
            % estimation of initial fit values
            A = pks(i); 
            % t0FirstEst = locs(i)-10;
            % [~,t0FirstEstPx] = min( abs( t-t0FirstEst ) );
            t0FirstEstPx = find(~prof_t_evnts_m(1:locPeaks_s(idx_p)), 1, 'last') - ...
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
                t0_est = locs(i)-10;
                tR_est = 3;
            end
            switch fitFun
                case 'expRise'
                    % parameters [t0, tR, A, y0]
                    x0 = [t0_est tR_est A-y0 y0];
                case 'sigmoid'
                    % parameters [bs A m tau]
                    x0 = [y0 A t0_est+(locs(i)-t0_est)/2 tR_est];
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












