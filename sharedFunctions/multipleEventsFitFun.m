function out = multipleEventsFitFun(t,y,ys,N_events,pos_peak,val_peak,baselineModel,...
    percBaseline,fitFun,t_ups,fitAproach,signOfPeak,mEvents,mBaseline,ax,useWeights)

% parameters of function:
% t = time
% y = data to fit
% ys = smoothed data of y
% N_events = number of peaks 
% pos_peak = position of peaks (points, not in time)
% val_peak = amplitude of peaks
% baselineModel = fit function for baseline (poly1, poly2, exp1, ... matlab build-in models)
% percBaseline = how many percent of signal is baseline
% fitFun = fit function to fit individual events
% fitAproach:    1 - fit baseline and all peaks together
%                2 - fit baseline and after fit peaks separatelly
%                3 - refit baseline and all peaks after fitted separatelly
% signOfPeak = positive or negative peak ('+' ,'-')
% mEvents = mask of events
% mBaseline = mask of baseline
% ax = handle to axes where the fit is plotted
% useWeights = use weighted fit 0 or 1

% make column vectors
t = t(:);
y = y(:);
ys = ys(:);

%keyboard

% set up positions and amplitudes of peaks    
[pos_peak,val_peak,~,~,~,~,~,~] = ...
    setUpEventFitFuntion(t,ys,N_events,pos_peak,val_peak,fitFun,signOfPeak);

% set up masks of events and baseline
[mEvents,mBaseline] = maskOfEvents(t,ys,N_events,pos_peak,mEvents,mBaseline,percBaseline);
 
% fit baseline with polynom
[~,fit_baseline,~,~] = fitOfBaseline(t,y,mBaseline,baselineModel,0);

% fit baseline with polynom
[~,fit_baselineS,~,~] = fitOfBaseline(t,ys,mBaseline,baselineModel,0);

% normalize to deltaF/F0, based on baseline fit
y = (y ./ fit_baseline) - 1;
ys = (ys ./ fit_baselineS) - 1;

% plot normalized profile
delete(findall(ax,'Tag','fitProfile'))
line(t,y,...
    'Parent',ax,'Color','k','LineStyle','-','LineWidth',1,...
    'Tag','fitProfile');
% set axes limits
ax.XLim = [min(t) max(t)];
ax.YLim = [min(y)*0.95 max(y)*1.05];

% set up positions and amplitudes of peaks, get initial guesses for fitting    
[pos_peak,val_peak,fun_event,p0,lb,ub,n_coeff_event,coeff_n] = ...
    setUpEventFitFuntion(t,ys,N_events,[],[],fitFun,signOfPeak);

% fit baseline with polynom
[f_baseline,fit_baseline,funBaseline,nParamsBaseline] = fitOfBaseline(t,y,mBaseline,baselineModel,1);


% figure
% plot(t,y,'ok')
% hold on 
% plot(t(mEvents),y(mEvents),'or')
% hold on 
% plot(t(mBaseline),y(mBaseline),'og')
% hold on 
% plot(t,fit_baseline,'m')



% do fitting
switch fitAproach
    
    case 1
        % fit baseline and all peaks together
        
        mEvents = true(size(mBaseline));
     
        events_coeff = fitOfData(p0,t,y,n_coeff_event,fun_event,N_events,...
            fit_baseline,f_baseline,mEvents,signOfPeak,fitFun,ax,useWeights,lb,ub);
        % delete partial fits
        delete(findobj(ax,'Tag','finalPartialFit'))
        p0 = [events_coeff,coeffvalues(f_baseline)];
        allCoeff = fitOfData(p0,t,y,n_coeff_event,fun_event,N_events,...
            fit_baseline,f_baseline,mEvents,signOfPeak,fitFun,ax,useWeights,...
            lb,ub,funBaseline,nParamsBaseline);
        
        events_coeff = allCoeff(1:end-nParamsBaseline);
        baseline_coeff = allCoeff(end-nParamsBaseline+1:end);
        fit_baseline = funBaseline(baseline_coeff,t);
        
    case 2
        % fit baseline and peaks separatelly
        events_coeff = fitOfData(p0,t,y,n_coeff_event,fun_event,N_events,...
            fit_baseline,f_baseline,mEvents,signOfPeak,fitFun,ax,useWeights,lb,ub);
        
        baseline_coeff = coeffvalues(f_baseline);
               
    case 3
        mEvents = true(size(mBaseline));
        keyboard
        % refit baseline and all peaks after fitted separatelly
        events_coeff = fitOfData(p0,t,y,n_coeff_event,fun_event,N_events,...
            fit_baseline,[],mEvents,signOfPeak,fitFun,ax,useWeights,lb,ub);

        p0 = [events_coeff,coeffvalues(f_baseline)];
        mEvents = true(size(mBaseline));
        
        allCoeff = fitOfData(p0,t,y,n_coeff_event,fun_event,N_events,...
            fit_baseline,f_baseline,mEvents,signOfPeak,fitFun,ax,useWeights,lb,ub,funBaseline,nParamsBaseline);
        
        events_coeff = allCoeff(1:end-nParamsBaseline);
        baseline_coeff = allCoeff(end-nParamsBaseline+1:end);
        fit_baseline = funBaseline(baseline_coeff,t);
end

% calculate final fit and fits of events
[F,F_events,F_individualEvents] = sumEvents(events_coeff,t,n_coeff_event,fun_event,N_events,fit_baseline);

% delete partial fits
delete(findobj(ax,'Tag','finalPartialFit'))


[F,F_events,F_individualEvents] = sumEvents(p0,t,n_coeff_event,fun_event,N_events,fit_baseline);

% figure
% plot(t,y)
% hold on
% plot(t,F,'r')


% create final output structure
events_n = num2cell(zeros(1,N_events));
coefficientsOfEvents = zeros(N_events,n_coeff_event);
for i=1:N_events
    events_n(1,i) = {sprintf('event %d',i)};
    coefficientsOfEvents(i,:) = events_coeff( 1+(i-1)*n_coeff_event : n_coeff_event+(i-1)*n_coeff_event );
end

out.t = t;
out.wholeFit = F;
out.allEventsFit = F_events;
out.baselineFit = fit_baseline;
out.individualEventsFits = [ events_n ; num2cell(F_individualEvents) ];
out.coefficientsOfFittedEvents = [[{''},coeff_n];[events_n',num2cell(coefficientsOfEvents)]];
out.bs = median(fit_baseline);
out.t0 = [];
out.yNorm = y;
out.ysNorm = ys;

if exist('t_ups', 'var')
    t_ups = t_ups(:);
    % evaluate events with different time axes
    [F_ups,F_events_ups,F_individualEvents_ups] = sumEvents(events_coeff,t_ups,n_coeff_event,fun_event,N_events,funBaseline(baseline_coeff,t_ups));
    out.t_ups.wholeFit = F_ups;
    out.t_ups.allEventsFit = F_events_ups;
    out.t_ups.baselineFit = funBaseline(baseline_coeff,t_ups);
    out.t_ups.individualEventsFits = [ events_n ; num2cell(F_individualEvents_ups) ];
    out.t_ups.t_ups = t_ups;
end

%toc
end


%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS USED IN FITTING %%%%%%%%%%%%%%%%%%%%%%%%
% set up event fit function; get initial guesses for fitting
function [pos_peak,val_peak,fun_event,p0,lb,ub,n_coeff_event,coeff_n] = ...
    setUpEventFitFuntion(t,ys,N_events,pos_peak,val_peak,fitFun,signOfPeak)

dt = mean(diff(t));

% incorporate also peak sign, in p0, lb, ub also fit function
switch signOfPeak
    case '+'
        ys = ys;
    case '-'
        keyboard
        ys = ys;
        % subtract 90 percentile as baseline, 
        % flip profile around zero .*(-1)
        % detrend, use streched exponential func to fit baseline
end

% get info about peaks of events
if isempty(pos_peak)
    % no peaks selected previously
    [val_peak,pos_peak,w,p] = findpeaks(ys,'NPeaks',N_events,'SortStr','descend');
else    
    % find all peaks and then take only those which amplitudes 
    % are same as amplitudes of peaks detected before 
    [pks,locs,w,p] = findpeaks(ys,'SortStr','descend');
    for i=1:N_events
        [~,ind(i)] = min(abs(locs-pos_peak(i)));
    end
    
    w = w(ind);
    p = p(ind);  
end

% sort based on location
[pos_peak,sI] = sort(pos_peak,'ascend');
val_peak = val_peak(sI);
w = w(sI);
p = p(sI);

% get initial guess for fitting
t0 = (pos_peak(:) - round(w(:)./2)).*dt; % in time units, ms
tpeak = pos_peak(:).*dt; % in time units, ms
F0 = zeros(size(t0));
tauR = (w(:).*dt)./3; % in time units, ms
tauD = w(:).*dt; % in time units, ms
A = val_peak(:);

switch fitFun
    
    case 'expModGauss'
        % function from: chromatography journal
        % convolution of CDF of normal distribution and exponential decay function     
        % coefficients: [A,m,sd,tau]
        fun_event = @(p,t) p(1).*exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - t./p(4)).*cdf('Normal',t,p(2)+(p(3)^2/p(4)),p(3));
        
        n_coeff_event = 4;
        coeff_n = {'A','m','sd','tauD'};
        switch signOfPeak
            case '+'
                p0 = [A,tpeak,w(:)./2,tauD];
                p0 = reshape(p0',[1,N_events*n_coeff_event]);
                lb = repmat([0 0 0 0],[1,N_events]);
                ub = repmat([inf inf inf inf],[1,N_events]);
            case '-'
                keyboard
        end

    
    case 'doubleBoltzmann'
        % fit with double boltzmann function
        % coefficients: [A, xH1, k1, xH2, k2]
        fun_event = @(p,x) p(1) .* ( 1./(1+exp((x-p(2))./p(3))) ) .* ( 1./(1+exp((x-p(4))./p(5))) ) ;
        
        n_coeff_event = 5;
        coeff_n = {'A','xH1','k1','xH2','k2'};
        switch signOfPeak
            case '+'
                % find closest point in t
                for i=1:N_events
                    [~,xH1(i,1)]= min(abs(t - (tpeak(i)-w(i)/2)));
                    [~,xH2(i,1)]= min(abs(t - (tpeak(i)+w(i)/2)));
                end
                p0 = [A,t(xH1),w(:)./2,t(xH2),-w(:)./2];
                p0 = reshape(p0',[1,N_events*n_coeff_event]);
                lb = repmat([0 0 0 0 -inf],[1,N_events]);
                ub = repmat([inf inf inf inf 0],[1,N_events]);
            case '-'
               keyboard 
        end
        
        
    case 'Asym2Sigmoid'
        % coefficients:[A, xc, w1, w2, w3]; originLab
        % Meanings: xc = center, A = amplitude, w1 = full width of half maximum,
        % w2 = variance of low-energy side, w3 = variance of high-energy side.
        % Lower Bounds: w1 > 0, w2 > 0, w3 > 0
        fun_event = @(p,x) p(1) .* ( 1./(1+exp(-(x-p(2)+p(3)/2)./p(4))) ) .* ( 1 - 1./(1+exp(-(x-p(2)-p(3)/2)./p(5))) ) ;
         
        n_coeff_event = 5;
        coeff_n = {'A','xc','w1','w2','w3'};
        switch signOfPeak
            case '+'
                p0 = [A,tpeak,w(:),w(:)/2,w(:)/2];
                p0 = reshape(p0',[1,N_events*n_coeff_event]);
                lb = repmat([0 0 0 0 0],[1,N_events]);
                ub = repmat([inf inf inf inf inf],[1,N_events]);
            case '-'
                keyboard
        end
        
        
    case '1expR1expD'
        % coefficients: [t0, tmax, tauR, A, tauD]
        fun_event = @(p,x) (x<p(1)).*0 + ...
                           (x>=p(1) & x<p(2)).*( (1-exp(-(x-p(1))./p(3))).* p(4) ) + ...
                           (x>=p(2)).*( (1-exp(-(p(2)-p(1))./p(3))).*exp(-(x-p(2))./p(5)) ).* p(4);   
        
                       
                       
%                        
%          fun_event = @(p,t) (t>=p(1)).*((1-exp(-(t-p(1))./p(3))) .*p(4) .*(exp(-(t-p(2))./p(5)))) + ...
%                             (t<p(1)).*0;

        n_coeff_event = 5;
        coeff_n = {'t0','tmax','tauR','A','tauD'};             
        switch signOfPeak
            case '+'
                p0 = [t0,tpeak,tauR,A,tauD];
                p0 = reshape(p0',[1,N_events*n_coeff_event]);
                lb = repmat([0 0 0 0 0],[1,N_events]);
                ub = repmat([inf inf inf inf inf],[1,N_events]);
            case '-'
                keyboard
        end               
           
%                     keyboard                
% figure
% plot(fun_event(p0,t))

    case '1expR2expD'
%         keyboard
        % coefficients: [t0, tmax, tauR, A, tauD1, tauD2]                                      
%         fun_event = @(p,t) (t<p(1)).*0 + ...
%             (t>=p(1)).*((1-exp(-(t-p(1))./p(3))) .*p(4) .*(exp(-(t-p(2))./p(5))) .*(exp(-(t-p(2))./p(6))));

                       
        fun_event = @(p,x) (x<p(1)).*0 + ...
            (x>=p(1) & x<p(2)).*( (1-exp(-(x-p(1))./p(3))).* p(4) ) + ...
            (x>=p(2)).*( (1-exp(-(p(2)-p(1))./p(3))).*exp(-(x-p(2))./p(5)).*(exp(-(t-p(2))./p(6))) ).* p(4);
        
                       
        n_coeff_event = 6;
        coeff_n = {'t0','tmax','tauR','A','tauD1','tauD2'};             
        switch signOfPeak
            case '+'
                p0 = [t0,tpeak,tauR,A,tauD,2*tauD];
                p0 = reshape(p0',[1,N_events*n_coeff_event]);
                lb = repmat([0 0 0 0 0 0],[1,N_events]);
                ub = repmat([inf inf inf inf inf inf],[1,N_events]);
            case '-'
                keyboard
        end               
           
%                     keyboard                
%  figure
%  plot(fun_event(p0,t))    
        
                                         
    case 'CaSpikeFun'
        %coefficients: [t0, FM, tA, tI, FI]
        fun_event = @(p,t) (t < p(1)) .* 0 + ...
                           (t >= p(1)) .* (0 + (p(2)*((-3*p(3)*(2 + (-1 + p(5))*p(3))*(4*p(3) - 7*p(4))*(p(3) - 3*p(4)))./(2*exp((2*(t-p(1)))/p(3)))+...
                                (3*p(3)*(p(3) - 3*p(4))*(p(3) - 2*p(4)))./exp((t-p(1))/p(3)) +(3*p(3)*(5 + (-1 + p(5))*p(3))*(2*p(3) - 5*p(4))*(p(3) - p(4)))./(5*exp((5*(t-p(1)))/p(3)))-...
                                (p(3)*(6 + (-1 + p(5))*p(3))*(p(3) - 2*p(4))*(p(3) - p(4)))./(6*exp((6.*(t-p(1)))/p(3))) -(6*p(4)^3*(1 + (-1 + p(5))*p(4)))./exp((t-p(1))/p(4))+...
                                (18*p(4)^3*(p(3) + p(4) + (-1 + p(5))*p(3)*p(4)))./(exp((p(3)^(-1) + p(4)^(-1)).*(t-p(1)))*(p(3) + p(4))) -... 
                                (18*p(4)^3*(p(3) + (2 + (-1 + p(5))*p(3))*p(4)))./(exp((2/p(3) + p(4)^(-1)).*(t-p(1)))*(p(3) + 2*p(4)))+...
                                (6*p(4)^3*(p(3) + (3 + (-1 + p(5))*p(3))*p(4)))./(exp((3/p(3) + p(4)^(-1)).*(t-p(1)))* (p(3) + 3*p(4))) -... 
                                (3*p(3)*(4 + (-1 + p(5))*p(3))*(5*p(3)^2 - 20*p(3)*p(4) +17*p(4)^2))./(4.*exp((4.*(t-p(1)))/p(3))) +...
                                (p(3)*(p(3)^2*(57 + 10*(-1 + p(5))*p(3)) - 3*p(3)*(84 + 13*(-1 + p(5))*p(3))*p(4) + (249 + 29*(-1 + p(5))*p(3))*p(4)^2))./(3.*exp((3.*(t-p(1)))/p(3)))-... 
                                ((-1 + p(5))*(p(3) - 3*p(4))*(p(3) - 2*p(4))*(p(3) - p(4))*(37*p(3)^4 + 252*p(3)^3*p(4) + 605*p(3)^2*p(4)^2 + 660*p(3)*p(4)^3 + 360*p(4)^4))./...
                                (60*(p(3) + p(4))*(p(3) + 2*p(4))*(p(3) + 3*p(4))) +  (6*(-1 + p(5))*p(3)^2*(p(3) - 3*p(4))*(p(3) - 2*p(4)).*...
                                cosh((t-p(1))/p(3)))./exp((2.*(t-p(1)))/p(3))))/...
                                ((p(3) - 3*p(4))*(p(3) - 2*p(4))*(p(3) - p(4))));
%          

                            
        n_coeff_event = 5;
        coeff_n = {'t0','FM','tA','tI','FI'};    
        switch signOfPeak
            case '+'
                p0 = [t0,A,tauR,tauD,ones(size(t0))];
                p0 = reshape(p0',[1,N_events*n_coeff_event]);
                lb = repmat([0 0 0 0 1],[1,N_events]);
                ub = repmat([inf inf inf inf 1],[1,N_events]);
            case '-'
                keyboard
        end
        
end

end


% set up masks of events and baseline
function [mEvents,mBaseline] = maskOfEvents(t,y,N_events,pos_peak,mEvents,mBaseline,percBaseline)

% sort peaks of events based on position of peak
pos_peak = sort(pos_peak,'ascend');

if ~isempty(mEvents)
    mEvents = logical(mEvents);
    
    if ~isempty(mBaseline)
        mBaseline = logical(mBaseline);
    else
        mBaseline = ~mEvents;
    end
    
else
    if ~isempty(percBaseline)
        % there is define treshhold for baseline (percentile)
        perc = prctile(y,percBaseline);
        m_exc = y > perc;
        
        mEvents = zeros(size(y));
        for i=1:numel(pos_peak)
            if m_exc(pos_peak(i))
                
                ind_e = find(m_exc(pos_peak(i):end)==0,1,'first')+pos_peak(i)-2;
                ind_s = find(m_exc(1:pos_peak(i))==0,1,'last') + 1;
                mEvents(ind_s:ind_e,1) = ones(ind_e-ind_s+1,1)*5;
                
            end
        end
        
    else
        % find events using orientation of 1.derivative
        mEvents = zeros(size(y));
        mBaseline = ones(size(y));
        for i=1:N_events
            if i==1
                if i==N_events
                    yp = y;
                    nPointsBeforePeak = numel(y(1:pos_peak(i)-1));
                    nPointsAfterPeak = numel(y(pos_peak(i)+1:end));
                else
                    yp = [y(1:pos_peak(i+1)) ; zeros(size(y(pos_peak(i+1)+1:end)))];
                    nPointsBeforePeak = numel(y(1:pos_peak(i)-1));
                    nPointsAfterPeak = floor(numel(y(pos_peak(i)+1:pos_peak(i+1)-1))*0.66);
                end
                
            elseif i==N_events && i~=1
                yp = [zeros(size(y(1:pos_peak(i-1)-1))) ; y(pos_peak(i-1):end)];
                nPointsBeforePeak = floor(numel(y(pos_peak(i-1)+1:pos_peak(i)-1))*0.33);
                nPointsAfterPeak = numel(y(pos_peak(i)+1:end));
                
            else
                yp = [zeros(size(y(1:pos_peak(i-1)-1))) ; y(pos_peak(i-1):pos_peak(i+1)) ; ...
                    zeros(size(y(pos_peak(i+1)+1:end)))];
                nPointsBeforePeak = floor(numel(y(pos_peak(i-1)+1:pos_peak(i)-1))*0.33);
                nPointsAfterPeak = floor(numel(y(pos_peak(i)+1:pos_peak(i+1)-1))*0.66);
            end
            
            if isempty(nPointsBeforePeak) nPointsBeforePeak=0; end
            if isempty(nPointsAfterPeak) nPointsAfterPeak=0; end

            % calculate gradient
            g = gradient(yp);
            g1 = flipud(g(1:pos_peak(i)-1));
            g2 = g(pos_peak(i)+1:end);
            
            ind_01 = find(g1 < 0, 1, 'first');
            ind_02 = find(g2 > 0, 1, 'first');
            ind_01 = pos_peak(i) - ind_01;
            ind_02 = pos_peak(i) + ind_02;
            
            if ~isempty(ind_01)
                ind_s = ind_01 - 1;
            end
            
            if ~isempty(ind_02)
                ind_e = ind_02 - 1;
            end
            
            mBaseline(ind_s:ind_e,1) = zeros(ind_e-ind_s+1,1);
            
            % compare two sets of indices, and assign broader interval for
            % spark
            if (pos_peak(i)-nPointsBeforePeak) < ind_s
                ind_s = pos_peak(i)-nPointsBeforePeak;
                if ind_s<1
                    ind_s = 1;
                end
            end
            
            if (pos_peak(i)+nPointsAfterPeak) > ind_e
                ind_e = pos_peak(i)+nPointsAfterPeak;
                if ind_e > numel(y) 
                    ind_e = numel(y);
                end
            end
                        
            mEvents(ind_s:ind_e,1) = ones(ind_e-ind_s+1,1);
            
            clearvars yp g g1 g2

        end
        
        mEvents = logical(mEvents);
        mBaseline = logical(mBaseline);
    end

end

end


% fit baseline
function [f_baseline,fit_baseline,funBaseline,nParamsBaseline] = fitOfBaseline(t,y,mBaseline,baselineModel,getCoeff)

if ~isempty(baselineModel)
    % fit
    switch baselineModel
        case 'stretchedExp'
            keyboard
            fitFun = @(p,x) p(1).*exp( (-x./p(2)).^p(3) ) + p(4);
            fitFunSum = @(p,t,y) sum( (y-fitFun(p,t)).^2 );
            % do fit
            opt = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
                    'MaxIter',1000,'MaxFunEvals',3000);
            coef = fmincon(@(p)fitFunSum(p,t(mBaseline),y(mBaseline)),...
                [1 -100 0.5 median(y(mBaseline))],[],[],[],[],[-inf -inf 0 -inf],[inf inf 1 inf],[],opt);
            
           figure
           plot(t,y)
           hold on
           plot(t,fitFun(coef,t),'r')
            fit_baseline = feval(f_baseline,t);
            
        otherwise
            ft = fittype(baselineModel);
            opts = fitoptions('Method','LinearLeastSquares',...
                'Normalize','on',...
                'Robust','Bisquare',...
                'Exclude',~mBaseline);
            
            [f_baseline, ~,~] = fit(t, y, ft, opts);
            fit_baseline = feval(f_baseline,t);
            
            % weights for next fit, push it to baseline
            wB = zeros(size(y));
            wB(mBaseline) = y(mBaseline) - fit_baseline(mBaseline);
            wB(wB<0)=0;
            wB = wB./(max(wB));
            wB(wB<=0)=1;
            % refit with weights
            opts.Weights = wB(mBaseline);
            opts.Exclude = [];
            [f_baseline, ~,~] = fit(t(mBaseline), y(mBaseline), ft, opts);
            fit_baseline = feval(f_baseline,t);
            
            if getCoeff
                % refit to get real parammeters of fit
                opts.Normalize = 'off';
                opts.Weights = [];
                opts.Exclude = [];
                opts.Robust = 'off';
                [f_baseline, ~,~] = fit(t, fit_baseline, ft, opts);
                fit_baseline = feval(f_baseline,t);
            end
            
            % cc = ( coeffvalues(f_baseline) - mean(t(mBaseline)) ) ./
            % std(t(mBaseline));
%             keyboard
%             
%             figure
%             plot(t(mBaseline),y(mBaseline),'ok')
%             hold on 
%             plot(t(mBaseline),wB(mBaseline),'or')
%             hold on 
%             plot(t,fit_baseline,'g')

            
            
    end
    
    % save baseline function
    switch baselineModel
        case 'exp1'
            funBaseline = @(p,x) p(1).*exp(p(2).*x);
            nParamsBaseline = 2;
        case 'exp2'
            funBaseline = @(p,x) p(1).*exp(p(2).*x) + p(3).*exp(p(4).*x);
            nParamsBaseline = 4;  
        case 'stretchedExp'
            funBaseline = @(p,x) p(1).*exp( (-x./p(2)).^p(3) ) + p(4);
            nParamsBaseline = 3;
        case 'poly1'
            funBaseline = @(p,x) p(1).*x + p(2);
            nParamsBaseline = 2;
        case 'poly2'
            funBaseline = @(p,x) p(1).*x.^2 + p(2).*x + p(3);
            nParamsBaseline = 3;
        case 'poly3'
            funBaseline = @(p,x) p(1).*x.^3 + p(2).*x.^2 + p(3).*x + p(4);
            nParamsBaseline = 4;
        case 'poly4'
            funBaseline = @(p,x) p(1).*x.^4 + p(2).*x.^3 + p(3).*x.^2 + p(4).*x + p(5);
            nParamsBaseline = 5;
        case 'poly5'
            funBaseline = @(p,x) p(1).*x.^5 + p(2).*x.^4 + p(3).*x.^3 + p(4).*x.^2 + p(5).*x + p(6);
            nParamsBaseline = 6;
        case 'poly6'
            funBaseline = @(p,x) p(1).*x.^6 + p(2).*x.^5 + p(3).*x.^4 + p(4).*x.^3 + p(5).*x.^2 + p(6).*x + p(7);
            nParamsBaseline = 7;
        case 'poly7'
            funBaseline = @(p,x) p(1).*x.^7 + p(2).*x.^6 + ...
                    p(3).*x.^5 + p(4).*x.^4 + p(5).*x.^3 + p(6).*x.^2 + p(7).*x + p(8);
            nParamsBaseline = 8;
        case 'poly8'
            funBaseline = @(p,x) p(1).*x.^8 + p(2).*x.^7 + p(3).*x.^6 + ...
                    p(4).*x.^5 + p(5).*x.^4 + p(6).*x.^3 + p(7).*x.^2 + p(8).*x + p(9);
            nParamsBaseline = 9;
        case 'poly9'
            funBaseline = @(p,x) p(1).*x.^9 + p(2).*x.^8 + p(3).*x.^7 + p(4).*x.^6 + ...
                    p(5).*x.^5 + p(6).*x.^4 + p(7).*x.^3 + p(8).*x.^2 + p(9).*x + p(10);
            nParamsBaseline = 10;
                
    end
     
else
    fit_baseline = zeros(size(y));
    f_baseline = @(t) zeros(size(t));
end

end


% do fitting
function events_coeff = fitOfData(p0,t,y,n_coeff_event,fun_event,N_events,...
    fit_baseline,f_baseline,mEvents,signOfPeak,fitFun,ax,useWeights,lb,ub,funBaseline,nParamsBaseline)
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
    
    try
        e(a,1) = find(mEvents(s(a):end) == 0,1,'first') + s(a) - 2;
    catch
        e(a,1) = numel(mEvents);
    end
    
    ind = e(a,1) + 1;    
    a = a+1; 
end
% keyboard
% fit individual masked events
events_coeff = zeros(size(p0));
for i=1:numel(s)
    
    if i==1 && i==numel(s)
        if ~exist('funBaseline','var')
            p0_mask = true(size(p0));
            N_events_part = sum(p0_mask)/n_coeff_event;
        else
            p0_mask = true(size(p0));
            N_events_part = (sum(p0_mask)-nParamsBaseline)/n_coeff_event;
        end
        t_mask = mEvents;
        
    else
        indx = ismember(p0,t(s(i):e(i)));
        indx = reshape(indx,[n_coeff_event,N_events]);
        
        for k=1:N_events
            if sum(indx(:,k))>=1
                indx(:,k) = ones(size(indx(:,k)));
            end
        end
        
        p0_mask = reshape(indx,[1,n_coeff_event*N_events]);
        t_mask = mEvents;
        t_mask(1:s(i))=0;
        t_mask(e(i):end)=0; 
        N_events_part = sum(p0_mask)/n_coeff_event;
    end
    
    if N_events_part<1
        continue
    end
    
    % part of profile
    yp = y(t_mask);
    tp = t(t_mask);
    bp = fit_baseline(t_mask);
    
%     figure
%     plot(tp,yp)
              
    % calculate weights for part of profile, further form 1 the bigger
    % weights, adjust for negative peaks
    switch signOfPeak
        case '+'
            if useWeights
                wE = (yp+1).^2;
            else
                wE = ones(size(yp));
            end
            
        case '-'
            keyboard
            if useWeights
                wE = ones(size(yp));
            else
                wE = ones(size(yp));
            end  
    end
    
    %%%%% do fitting
    options = optimoptions('fmincon','TolFun',1e-12,'TolX',1e-12,'TolCon',1e-12,...
        'MaxIter',1000,'MaxFunEvals',3000,'UseParallel',0,'Display','none',...
        'OutputFcn',@(x,optimValues,state)outFcnFmincon(x,optimValues,state,ax,tp,n_coeff_event,fun_event,N_events_part,bp));
   
    if ~exist('funBaseline','var')

        events_coeff(p0_mask) = fmincon(@(p)fitFunSumOfSquaersW(p,tp,yp,n_coeff_event,fun_event,N_events_part,bp,wE),...
            p0(p0_mask),[],[],[],[],lb(p0_mask),ub(p0_mask),[],options);
    else
     % keyboard
        delete(findobj(ax,'Tag','finalPartialFit'))
        
        options = optimoptions('patternsearch','TolFun',1e-12,'TolX',1e-12,'TolCon',1e-12,'MeshTolerance',1e-12,...
        'MaxIter',1000,'MaxFunEvals',3000,'UseParallel',0,'Display','none',...
        'OutputFcn',@(x,optimValues,state)outFcnGlobalOpt(x,optimValues,state,ax,tp,n_coeff_event,fun_event,N_events_part,bp,funBaseline,nParamsBaseline));
    
        events_coeff(p0_mask) = ...
            patternsearch(@(p)fitFunSumOfSquaersW(p,tp,yp,n_coeff_event,fun_event,N_events_part,bp,wE,funBaseline,nParamsBaseline),...
            p0(p0_mask),[],[],[],[],...
            [lb(p0_mask(1:end-nParamsBaseline)),-inf([1,nParamsBaseline])],...
            [ub(p0_mask(1:end-nParamsBaseline)),inf([1,nParamsBaseline])],[],options);
           
        
%         gs = GlobalSearch;
% ff = @(p)fitFunSumOfSquaersW(p,tp,yp,n_coeff_event,fun_event,N_events_part,bp,wE,funBaseline,nParamsBaseline);
% problem = createOptimProblem('fmincon','x0',p0,...
%     'objective',ff,'lb',[lb(p0_mask(1:end-nParamsBaseline)),-inf([1,nParamsBaseline])],'ub',[ub(p0_mask(1:end-nParamsBaseline)),inf([1,nParamsBaseline])]);
% events_coeff = run(gs,problem);
        
        
        
    end
            
% % % @(p)nonlinEq(p,tp,fun_event,fitFun,n_coef_event,N_events_part)            
% % %         
% 
% gs = GlobalSearch;
% ff = @(p)fitFunSumOfSquaersW(p,tp,yp,n_coeff_event,fun_event,N_events_part,bp,wE,funBaseline,nParamsBaseline);
% problem = createOptimProblem('fmincon','x0',p0,...
%     'objective',ff,'lb',[lb(p0_mask(1:end-nParamsBaseline)),-inf([1,nParamsBaseline])],'ub',[ub(p0_mask(1:end-nParamsBaseline)),inf([1,nParamsBaseline])]);
% x = run(gs,problem);

% 
%    [F,~,~] = sumEvents(events_coeff,tp,n_coeff_event,fun_event,N_events_part,bp,funBaseline,nParamsBaseline);
% 
%              figure
%              plot(tp,yp)
%              hold on 
%              plot(tp,F,'r')
            
end
    
end


% sum N_events together using chosen function and proper initial
% parameters
function [F,F_events,F_individualEvents] = sumEvents(p,t,n_coef_event,fun_event,N_events,fit_baseline,funBaseline,nParamsBaseline)

F_events = zeros(numel(t),1);
F_individualEvents = zeros(numel(t),N_events);
for j=1:N_events
    
    F_event = fun_event(p(1+(j-1)*n_coef_event:n_coef_event+(j-1)*n_coef_event),t);
    F_event(~isfinite(F_event))=0;
    
    F_events = F_events + F_event;
    F_individualEvents(:,j) = F_event;
    
end

if ~exist('funBaseline','var')
    F = F_events + fit_baseline;
else
    F = F_events + funBaseline(p(end-nParamsBaseline+1:end),t);
end
    
end

% fit function to use with fminsearch, fmincon; function returning
% scalar % weighted
function resSum = fitFunSumOfSquaersW(p,t,y,n_coeff_event,fun_event,N_events,fit_baseline,wE,funBaseline,nParamsBaseline)

if ~exist('funBaseline','var')
    [F,~,~] = sumEvents(p,t,n_coeff_event,fun_event,N_events,fit_baseline);
else
    [F,~,~] = sumEvents(p,t,n_coeff_event,fun_event,N_events,fit_baseline,funBaseline,nParamsBaseline);
end

resSum = sum(wE.*(y - F).^2);

end

% nonlinear equalities for fmincon
function [c,ceq] = nonlinEq(p,x,fun,model,n_coef_event,N_events_part)
%keyboard
% delta x
dt = mean(diff(x));

% break points and equalities
switch model
           
    case 'expModGaussPiecewise'
         %keyboard            
        % equal value of function in break point 
        c = [];
        ceq = [];
        
    case '1expR1expD'
                          
        br = p(1:n_coef_event:N_events_part*n_coef_event);
               
        % equal value of function in break point 
        c = [];
        ceq(1:numel(br)) = fun(p,br-dt) - fun(p,br);
        ceq = [];
                 

    otherwise        
        % equal F
        c = [];
        ceq = [];
                       
end

% % equal value of differential of function in break point
%ceq(2) = diff(fun(p,(br:dt:br+dt))) - diff(fun(p,(br-dt:dt:br)));

end


% output function for fmincon
function stop = outFcnFmincon(x,optimValues,state,ax,tp,n_coef_event,fun_event,N_events_part,bp)

stop = false;

% stop fitting if stop button in fit profile window is pressed
stop = getappdata(ax.Parent,'stopFiting');

[F,~,~] = sumEvents(x,tp,n_coef_event,fun_event,N_events_part,bp);
      
%delete(findobj(ax,'Tag','partialFit'))

hl = findobj(ax,'Tag','partialFit');

switch state
    
    case 'iter'
        %keyboard
        % Make updates to plot or guis as needed
        set(hl,'XData',tp,'YData',F)
%         line(tp,F,'Parent',ax,'Color','b','LineStyle','-','LineWidth',2,...
%             'Tag','partialFit');
        drawnow
        
    case 'interrupt'
        % Probably no action here. Check conditions to see
        % whether optimization should quit.
        
    case 'init'
        % Setup for plots or guis             
        line(tp,F,'Parent',ax,'Color','b','LineStyle','-','LineWidth',2,...
                'Tag','partialFit');
        drawnow
        
    case 'done'
        % Cleanup of plots, guis, or final plot  
        set(hl,'XData',tp,'YData',F,'Tag','finalPartialFit')
        
%         line(tp,F,'Parent',ax,'Color','b','LineStyle','-','LineWidth',2,...
%             'Tag','finalPartialFit')
        drawnow
        
    otherwise
        
end

end


% output function for paternsearch
function [stop,options,optchanged] = outFcnGlobalOpt(x,optimValues,state,ax,tp,n_coef_event,fun_event,N_events_part,bp,funBaseline,nParamsBaseline)

stop = false;

% stop fitting if stop button in fit profile window is pressed
stop = getappdata(ax.Parent,'stopFiting');

options = [];
optchanged = [];

[F,~,~] = sumEvents(x.x,tp,n_coef_event,fun_event,N_events_part,bp,funBaseline,nParamsBaseline);
      
%delete(findobj(ax,'Tag','partialFit'))

hl = findobj(ax,'Tag','partialFit');

switch state
    
    case 'iter'
        %keyboard
        % Make updates to plot or guis as needed
        set(hl,'XData',tp,'YData',F)
%         line(tp,F,'Parent',ax,'Color','b','LineStyle','-','LineWidth',2,...
%             'Tag','partialFit');
        drawnow
        
    case 'interrupt'
        % Probably no action here. Check conditions to see
        % whether optimization should quit.
        
    case 'init'
        % Setup for plots or guis             
        line(tp,F,'Parent',ax,'Color','b','LineStyle','-','LineWidth',2,...
                'Tag','partialFit');
        drawnow
        
    case 'done'
        % Cleanup of plots, guis, or final plot  
        set(hl,'XData',tp,'YData',F,'Tag','finalPartialFit')
        
%         line(tp,F,'Parent',ax,'Color','b','LineStyle','-','LineWidth',2,...
%             'Tag','finalPartialFit')
        drawnow
        
    otherwise
        
end

end












