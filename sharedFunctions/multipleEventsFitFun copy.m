function out = multipleEventsFitFun(pos_peak,val_peak,t,y,baselineModel,...
    percBaseline,fitFun,fitting,t_ups,recalculate,signOfPeak,mEvents,mBaseline)

%tic
% parameters of function
% pos_peak = position of peaks (not in time)
% val_peak = amplitude of peaks
% t = time
% y = data to fit
% baselineModel = fit function for baseline (poly1, poly2 ...; matlab build in models)
% percBaseline = how many percent of signal is baseline
% fitFun = fit function to fit individual events
% fitting = fminsearch, lsqnonlin, fmincon
% recalculate = refit baseline and peaks 0 or 1
% signOfPeak = positive or negative peak ('+' ,'-')
% mEvents = mask of events
% mBaseline = mask of baseline
keyboard
% make column vectors
t = t(:);
y = y(:);
N_events = numel(pos_peak);

if ~isempty(mEvents)
    
    mEvents = logical(mEvents);
    if ~isempty(mBaseline)
    
    mBaseline = logical(mBaseline);
    else
        mBaseline = ~mEvents;
    end
else
    
    perc = prctile(y,percBaseline);
    m_exc = y > perc;
    
    %%%%% fit only selected peaks from signal, peaks with known position
    % set weights of signal
    mEvents = zeros(size(y));
    for i=1:numel(pos_peak)
        if m_exc(pos_peak(i))
            
            ind_e = find(m_exc(pos_peak(i):end)==0,1,'first')+pos_peak(i)-2;
            ind_s = find(m_exc(1:pos_peak(i))==0,1,'last') + 1;
            %         ind_s_o = ind_s;
            %         p_peak = pos_peak(i) - ind_s + 2;
            %
            %         y_p = smooth(y(ind_s:ind_e),10/length(y(ind_s:ind_e)),'loess');
            %
            %         g = gradient(y_p);
            %         g1 = flipud(g(1:p_peak-1));
            %         g2 = g(p_peak+1:end);
            %
            % %
            % %         % keyboard
            % %         figure
            % %         plot(y(ind_s:ind_e))
            % %         hold on
            % %         plot(y_p,'r')
            % %         hold on
            % %         plot(g,'g')
            % %         hold on
            % %         plot(zeros(size(g)))
            %
            %         ind_01 = find(g1 < 0, 1, 'first');
            %         ind_02 = find(g2 > 0, 1, 'first');
            %         ind_01 = p_peak - ind_01;
            %         ind_02 = p_peak + ind_02;
            %
            %         if ~isempty(ind_01)
            %             ind_s = ind_s_o + ind_01 - 1;
            %         end
            %
            %         if ~isempty(ind_02)
            %             ind_e = ind_s_o + ind_02 - 1;
            %         end
            
            
            
            mEvents(ind_s:ind_e,1) = ones(ind_e-ind_s+1,1)*5;
            
            clearvars y_p g g1 g2
            
        end
    end
    
    mEvents = logical(mEvents);
    
end


%%%%% find values of negative peaks in SR signal, position of peaks taken from
% cytosolic signal
pos_peak = pos_peak(:);
switch signOfPeak
    
    case '+'
        val_peak = val_peak(:);
        
    case '-'
        
        for i=1:numel(pos_peak)
            % find minimum in SR signal around (+- 10 ms) peak position in
            % cytosolic signal
            [v,p] = min(y(pos_peak(i)-round(10/(t(2)-t(1))):pos_peak(i)+round(10/(t(2)-t(1)))));
            pos_peak(i,1) = pos_peak(i)-round(10/(t(2)-t(1)))+p-1;
            val_peak(i,1) = -v;
        end
        
end


% keyboard
% % najskor fitovat lsqnonlin kvoli rychlosti a az potom to fitnut fminconom
% % ak to padne, pridat odhad baseliny, nejaky sofistikovanejsi, detrend, PID
% % alebo nejaky filter, ARMA, prediction model ???
%
% %%%%%%%%%%%%%


% remove baseline artefact, fit baseline with polynom or smooth
[f_baseline,fit_baseline] = fitOfBaseline(t,y,mBaseline,baselineModel);

% event fit function
[fun_event,p0,n_coef_event,coeff_n] = setFitFunction(fitFun,N_events,t,val_peak,pos_peak);

% do fitting
events_coeff = fitOfData(p0,t,y,n_coef_event,fun_event,N_events,fit_baseline,fitting,mEvents,signOfPeak);

% calculate final fit and fits of events 
[F,F_events,F_individualEvents] = sumEvents(events_coeff,t,n_coef_event,fun_event,N_events,fit_baseline);

% recalculate baseline with removed fits and repeat fitting procedure
if recalculate
      
    % set weights, to fit properly peaks, bigger number = bigger weight       
    % mask of events   
       
    switch signOfPeak
        
        case '+'
            
            mEvents1 = mEvents(:) + F_events(:)>0.001;
            % set weights for events.
            wE = zeros(size(y));
            wE(~mEvents1) = 0.001;
            wE(mEvents1) = (y(mEvents1)+1).^2;
            
        case '-'
            mEvents1 = mEvents(:) + F_events(:)<-0.001;
            % set weights for events
            wE = zeros(size(y));
            wE(~mEvents1) = 1;
            wE(mEvents1) = 1;
            %wE(mEvents1) = (abs(y(mEvents1)-1)+1).^2;

    end    
        
    
%     keyboard
%     
%     figure
%     plot(y)
%     hold on 
%     plot(F,'r')
    
    
    [f_baseline,fit_baseline] = fitOfBaseline(t,y,mEvents1,baselineModel);
    
    % do fitting of whole profile 
    options = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
        'MaxIter',1000,'MaxFunEvals',3000,'UseParallel',0);
            
    events_coeff = fmincon(@(p)fitFunSumOfSquaersW(p,t,y,n_coef_event,fun_event,N_events,fit_baseline,wE),...
        events_coeff,[],[],[],[],zeros(1,numel(p0)),inf(1,numel(p0)),@(p)nonlinEq(p,x,fitFun,model),options);
    
    % calculate final fit and fits of events 
    [F,F_events,F_individualEvents] = sumEvents(events_coeff,t,n_coef_event,fun_event,N_events,fit_baseline);
    
end

% create final output structure
events_n = num2cell(zeros(1,N_events));
coefficientsOfEvents = zeros(N_events,n_coef_event);
for i=1:N_events
    events_n(1,i) = {sprintf('event %d',i)};
    coefficientsOfEvents(i,:) = events_coeff( 1+(i-1)*n_coef_event : n_coef_event+(i-1)*n_coef_event );
end

out.t = t;
out.wholeFit = F;
out.allEventsFit = F_events;
out.baselineFit = fit_baseline;
out.individualEventsFits = [ events_n ; num2cell(F_individualEvents) ];
out.coefficientsOfFittedEvents = [[{''},coeff_n];[events_n',num2cell(coefficientsOfEvents)]];

if exist('t_ups', 'var')
    t_ups = t_ups(:);
    % evaluate events with different time axes
    [F_ups,F_events_ups,F_individualEvents_ups] = sumEvents(events_coeff,t_ups,n_coef_event,fun_event,N_events,feval(f_baseline,t_ups));
    out.t_ups.wholeFit = F_ups;
    out.t_ups.allEventsFit = F_events_ups;
    out.t_ups.baselineFit = feval(f_baseline,t_ups);
    out.t_ups.individualEventsFits = [ events_n ; num2cell(F_individualEvents_ups) ];
    out.t_ups.t_ups = t_ups;
end

%toc
end


%%%%%%%%%%%%%%%%%%%%%%%% functions for fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit baseline
function [f_baseline,fit_baseline] = fitOfBaseline(t,y,mBaseline,baselineModel)

switch baselineModel
    
    case 'smooth'
        s_span = numel(y)/(numel(y)+ceil(5/(t(2)-t(1))));
        ft = fittype('smoothingspline');
        opts = fitoptions('Method','SmoothingSpline',...
                          'SmoothingParam',s_span,...
                          'Normalize','on',...
                          'Exclude',~mBaseline);
        % fit
        [f_baseline, ~,~] = fit(t, y, ft, opts);
        fit_baseline = feval(f_baseline,t);
    
    case 'fitSpline'
        keyboard
        
        
        figure
        plot(t(mBaseline),y(mBaseline))
        hold on 
        plot(t,fit_baseline,'r')
        
        % fit spline to profile, so far the best
        %splinetool(t_all,F_all)
        n_knots = ceil(t(end)/200);
        splOrd = 3; % spline order, poly3
        
        spl = spap2(n_knots,splOrd,t(~mEvents),y(~mEvents));
        % newknt for a possibly better knot distribution
        knots = newknt(spl);
%         [~,indxK] = min(abs(knots-startH2O2));
%         knots(indxK) = startH2O2;
        % least-squares approximation
        spl = spap2(knots, splOrd, t(~mEvents), y(~mEvents));
        
        fit_baseline = fnval(spl,t);
        % calculate 1.derivative of fit
        dF_all_fit = fnval(fnder(spl,1),t_all);
        %dF_all_fit = gradient(F_all_fit);
        
    otherwise      
        ft = fittype(baselineModel);
        opts = fitoptions('Method','LinearLeastSquares',...
                          'Normalize','on',...
                          'Robust','Bisquare',...
                          'Exclude',~mBaseline);
        % fit
        [f_baseline, ~,~] = fit(t, y, ft, opts);
        fit_baseline = feval(f_baseline,t);
end



end


% do fitting
function events_coeff = fitOfData(p0,t,y,n_coef_event,fun_event,N_events,fit_baseline,fitting,mEvents,signOfPeak)
% find starts and ends of parts of profile with events, they can be
% merged or just single events
a = 1;
ind = 1;
while ind < numel(t)
   
    try
        s(a,1) = find(mEvents(ind:end) == 1,1,'first') + ind - 1;
    catch
       break 
    end
    e(a,1) = find(mEvents(s(a):end) == 0,1,'first') + s(a) - 2;
    
    ind = e(a,1) + 1;    
    a = a+1; 
end

% fit individual masked events
events_coeff = zeros(size(p0));
for i=1:numel(s)
    
    indx = ismember(p0,t(s(i):e(i)));
    indx = reshape(indx,[n_coef_event,N_events]);
    
    for k=1:N_events
        if sum(indx(:,k))>=1
            
            indx(:,k) = ones(size(indx(:,k)));
        end
    end
    
    p0_mask = reshape(indx,[1,n_coef_event*N_events]);
    t_mask = mEvents;
    t_mask(1:s(i))=0;
    t_mask(e(i):end)=0;
    
    N_events_part = sum(p0_mask)/n_coef_event;
          
    % part of profile
    yp = y(t_mask);
    tp = t(t_mask);
    bp = fit_baseline(t_mask);
    
    % calculate weights for part of profile, further form 1 the bigger weight
    switch signOfPeak
        
        case '+'
            wE = (yp+1).^2;
            
        case '-'
            % keyboard
            wE = (yp-1).^2;
            
    end
    
    %wE = ones(size(yp));
    
    
    %%%%% do fitting
    switch fitting
        
        case 'simplex'
            try
                options = optimset('MaxFunEvals',30000,'MaxIter',10000,'TolFun',1e-12,'TolX',1e-12);
                events_coeff(p0_mask) = fminsearch(@(p)fitFunSumOfSquaersW(p,tp,yp,n_coef_event,fun_event,N_events_part,bp,wE),...
                                            p0(p0_mask),options);               
            catch               
                options = optimoptions('fmincon','TolFun',1e-12,'TolX',1e-12,'TolCon',1e-12,...
                    'MaxIter',10000,'MaxFunEvals',30000,'UseParallel',0);
                
                events_coeff(p0_mask) = fmincon(@(p)fitFunSumOfSquaersW(p,tp,yp,n_coef_event,fun_event,N_events_part,bp,wE),...
                    p0(p0_mask),[],[],[],[],zeros(1,numel(p0(p0_mask))),inf(1,numel(p0(p0_mask))),...
                    @(p)nonlinEq(p,tp,fun_event,'2expR2expD',n_coef_event,N_events_part),options);
                
            end
                    
        case 'fmincon'
                                 
            options = optimoptions('fmincon','TolFun',1e-12,'TolX',1e-12,'TolCon',1e-12,...
                'MaxIter',10000,'MaxFunEvals',30000,'UseParallel',0);
            
            events_coeff(p0_mask) = fmincon(@(p)fitFunSumOfSquaersW(p,tp,yp,n_coef_event,fun_event,N_events_part,bp,wE),...
                p0(p0_mask),[],[],[],[],zeros(1,numel(p0(p0_mask))),inf(1,numel(p0(p0_mask))),...
                @(p)nonlinEq(p,tp,fun_event,'2expR2expD',n_coef_event,N_events_part),options);
            
            %@(p)nonlinEq(p,tp,fun_event,'2expR2expD')
            
            
            
            
%             
%             [F,F_events,F_individualEvents] = sumEvents(events_coeff,tp,n_coef_event,fun_event,N_events,bp);
%             
%             
%             figure
%             plot(tp,yp)
%             hold on 
%             plot(tp,F,'r')
            
    end
    
end

end


% set fit function
function [fun_event,p0,n_coef_event,coeff_n] = setFitFunction(fitFun,N_events,t,val_peak,pos_peak)
switch fitFun
    
    case 'expModGauss'
        % function from:
        % Lacouture, Y., & Cousineau, D. (2008).
        % How to use MATLAB to fit the ex-Gaussian and other probability functions to a distribution of response times.
        % Tutorials in Quantitative Methods for Psychology, 4(1), 35?45. http://doi.org/10.20982/tqmp.04.1.p035
        % convolution of CDF of normal distribution and exponential decay function
        % coefficients: [t0,A,m,sd,tau]
        % fun_event = @(p,t) (t<p(1)).*0 + (t>=p(1)).*(p(2).*exp(p(3)/p(5) + p(4)^2/(2*p(5)^2) - t./p(5)).*cdf('Normal',t,p(3)+(p(4)^2/p(5)),p(4)));
        % coefficients: [A,m,sd,tau]
        fun_event = @(p,t) p(1).*exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - t./p(4)).*cdf('Normal',t,p(2)+(p(3)^2/p(4)),p(3));
        n_coef_event = 4;
        coeff_n = {'A','m','sd','tau'};
        
        p0 = zeros(1,N_events*n_coef_event);
        for i=1:N_events
            
            p0(1 + (i-1)*n_coef_event) = val_peak(i);
            p0(2 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(3 + (i-1)*n_coef_event) = 10;
            p0(4 + (i-1)*n_coef_event) = 20;
            
        end
        
    case 'exp2ModGauss'
        % function from:
        % Lacouture, Y., & Cousineau, D. (2008).
        % How to use MATLAB to fit the ex-Gaussian and other probability functions to a distribution of response times.
        % Tutorials in Quantitative Methods for Psychology, 4(1), 35?45. http://doi.org/10.20982/tqmp.04.1.p035
        % convolution of CDF of normal distribution and exponential decay function
        % coefficients: [t0,A,m,sd,tau]
        % fun_event = @(p,t) (t<p(1)).*0 + (t>=p(1)).*(p(2).*exp(p(3)/p(5) + p(4)^2/(2*p(5)^2) - t./p(5)).*cdf('Normal',t,p(3)+(p(4)^2/p(5)),p(4)));
        % coefficients: [A,m,sd,tau1,tau2]
        fun_event = @(p,t) p(1).*exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - t./p(4)).*exp(p(2)/p(5) + p(3)^2/(2*p(5)^2) - t./p(5)).*cdf('Normal',t,p(2)+(p(3)^2/p(4)),p(3));
        n_coef_event = 5;
        coeff_n = {'A','m','sd','tau1','tau2'};
        
        p0 = zeros(1,N_events*n_coef_event);
        for i=1:N_events
            
            p0(1 + (i-1)*n_coef_event) = val_peak(i);
            p0(2 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(3 + (i-1)*n_coef_event) = 10;
            p0(4 + (i-1)*n_coef_event) = 20;
            p0(5 + (i-1)*n_coef_event) = 100;
            
        end
        
        
    case '1expR1expD'
        %coefficients: [t0, tmax, tauR1, Ampl, tauD1]
        fun_event = @(p,t) (t>=p(1)).*((1-exp(-(t-p(1))./p(3))) .*p(4) .*(exp(-(t-p(2))./p(5)))) + ...
            (t<=p(1)).*0;
        n_coef_event = 5;
        coeff_n = {'t0','tmax','tauR1','Ampl','tauD1'};
        
        p0 = zeros(1,N_events*n_coef_event);
        for i=1:N_events
            
            p0(1 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(2 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(3 + (i-1)*n_coef_event) = 5;
            p0(4 + (i-1)*n_coef_event) = val_peak(i);
            p0(5 + (i-1)*n_coef_event) = 10;
            
        end
        
    case '1expR2expD'
        %coefficients: [t0, tmax, tauR1, Ampl, tauD1, tauD2]
        fun_event = @(p,t) (t>=p(1)).*((1-exp(-(t-p(1))./p(3))) .*p(4) .*(exp(-(t-p(2))./p(5)) + exp(-(t-p(2))./p(6)))) + ...
            (t<p(1)).*0;
        n_coef_event = 6;
        coeff_n = {'t0','tmax','tauR1','Ampl','tauD1','tauD2'};
        
        p0 = zeros(1,N_events*n_coef_event);
        for i=1:N_events
            
            p0(1 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(2 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(3 + (i-1)*n_coef_event) = 5;
            p0(4 + (i-1)*n_coef_event) = val_peak(i);
            p0(5 + (i-1)*n_coef_event) = 10;
            p0(6 + (i-1)*n_coef_event) = 20;
            
        end
        
    case '2expR2expD'
        %coefficients: [t0, tmax, tauR1, tauR2, Ampl, tauD1, tauD2]
        fun_event = @(p,t) (t>=p(1)).*(((1-exp(-(t-p(1))./p(3))).*(1-exp(-(t-p(1))./p(4)))) .*p(5) .*(exp(-(t-p(2))./p(6)) + exp(-(t-p(2))./p(7)))) + ...
            (t< p(1)).*0;
        n_coef_event = 7;
        coeff_n = {'t0','tmax','tauR1','tauR2','Ampl','tauD1','tauD2'};
        
        p0 = zeros(1,N_events*n_coef_event);
        for i=1:N_events
            
            p0(1 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(2 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(3 + (i-1)*n_coef_event) = 3;
            p0(4 + (i-1)*n_coef_event) = 5;
            p0(5 + (i-1)*n_coef_event) = val_peak(i);
            p0(6 + (i-1)*n_coef_event) = 10;
            p0(7 + (i-1)*n_coef_event) = 20;
            
        end
        
        
    case 'sigR1expD'
        keyboard 
        %coefficients: [t0, tmax, tauR1, tauR2, Ampl, tauD1, tauD2]
        fun_event = @(p,t) (t>=p(1)).*(((1-exp(-(t-p(1))./p(3))).*(1-exp(-(t-p(1))./p(4)))) .*p(5) .*(exp(-(t-p(2))./p(6)) + exp(-(t-p(2))./p(7)))) + ...
            (t< p(1)).*0;
        n_coef_event = 7;
        coeff_n = {'t0','tmax','tauR1','tauR2','Ampl','tauD1','tauD2'};
        
        p0 = zeros(1,N_events*n_coef_event);
        for i=1:N_events
            
            p0(1 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(2 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(3 + (i-1)*n_coef_event) = 3;
            p0(4 + (i-1)*n_coef_event) = 5;
            p0(5 + (i-1)*n_coef_event) = val_peak(i);
            p0(6 + (i-1)*n_coef_event) = 10;
            p0(7 + (i-1)*n_coef_event) = 20;
            
        end
        
        
    case 'ArdosLikeFun'
        %coefficients: [t0, u, tauR, A1, tauD1, A2, tauD2]
        fun_event = @(p,t) (t>=p(1)).*( (1-exp(-(t-p(1))./p(3))) .* p(4).*((exp(-(t-p(1)-p(2))./p(5))) + p(6).*(exp(-(t-p(1)-p(2))./p(7)))) ) + ...
            (t<p(1)).*0;
        n_coef_event = 7;
        coeff_n = {'t0','u','tR','A1','tD1','A2','tD2'};
        
        p0 = zeros(1,N_events*n_coef_event);
        for i=1:N_events
            
            p0(1 + (i-1)*n_coef_event) = t(pos_peak(i));
            p0(2 + (i-1)*n_coef_event) = 1;
            p0(3 + (i-1)*n_coef_event) = 3;
            p0(4 + (i-1)*n_coef_event) = val_peak(i);
            p0(5 + (i-1)*n_coef_event) = 10;
            p0(6 + (i-1)*n_coef_event) = val_peak(i);
            p0(7 + (i-1)*n_coef_event) = 20;
            
        end
        %%%% try convolution of 3 exponential functions
        
end

end


% sum N_events together using chosen function and proper initial
% parameters
function [F,F_events,F_individualEvents] = sumEvents(p,t,n_coef_event,fun_event,N_events,fit_baseline)

F_events = zeros(numel(t),1);
F_individualEvents = zeros(numel(t),N_events);
for j=1:N_events
    
    F_event = fun_event(p(1+(j-1)*n_coef_event:n_coef_event+(j-1)*n_coef_event),t);
    F_event(~isfinite(F_event))=0;
    
    F_events = F_events + F_event;
    F_individualEvents(:,j) = F_event;
    
end

F = F_events + fit_baseline;

end


% fit function to use with lsqnonlin; function returning
% vector
function diff = fitFunDiff(p,t,y,n_coef_event,fun_event,N_events,fit_baseline)

[F,~,~] = sumEvents(p,t,n_coef_event,fun_event,N_events,fit_baseline);
diff = y - F;

end


% fit function to use with fminsearch, fmincon; function returning
% scalar
function resSum = fitFunSumOfSquaers(p,t,y,n_coef_event,fun_event,N_events,fit_baseline)

[F,~,~] = sumEvents(p,t,n_coef_event,fun_event,N_events,fit_baseline);
resSum = sum((y - F).^2);

end
% weighted
function resSum = fitFunSumOfSquaersW(p,t,y,n_coef_event,fun_event,N_events,fit_baseline,wE)

[F,~,~] = sumEvents(p,t,n_coef_event,fun_event,N_events,fit_baseline);
resSum = sum(wE.*(y - F).^2);

end


% nonlinear equalities for fmincon
function [c,ceq] = nonlinEq(p,x,fun,model,n_coef_event,N_events_part)
%keyboard
% delta x
dt = mean(diff(x));

% break points and equalities
switch model
        
    case '2expR2expD'
                  
        br = p(1:n_coef_event:N_events_part*n_coef_event);
               
        % equal value of function in break point 
        c = [];
        ceq(1:numel(br)) = fun(p,br-dt) - fun(p,br);
%         ceq(1) = fun(p,br2-dt) - fun(p,br2);
        
    otherwise        
        % equal F
        c = [];
        ceq = [];
                       
end

% % equal value of differential of function in break point
%ceq(2) = diff(fun(p,(br:dt:br+dt))) - diff(fun(p,(br-dt:dt:br)));

end


















