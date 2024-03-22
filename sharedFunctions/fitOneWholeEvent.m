function out = fitOneWholeEvent(x, y, model, options)
%{
parameters of function
x = independent variable
y = data to fit 
model:
    1. 1expR1expD (single exponential rise and decay)
    2. spline
    3. constThenSpline (constant for x<t0 and spline for x>=t0)
    4. EMG (exponentially modified gaussian)
    5. CaSpike model
    6. Gauss
    7. doubleBoltzmann
    8. Asym2Sigmoid
options:
    t0 = start of event
    bs = baseline before event
    x_ups = upscaled independent variable
    mProf = mask of profile
    splineOrder
    spline_knots_seq = sequence of knots
    addMoreWeightToEvent
%}
arguments
    x (:,1) {mustBeNumeric}            
    y (:,1) {mustBeNumeric}             
    model (1,1) {mustBeText, ...
        mustBeMember(model,{'1expR1expD','spline', ...
        'constThenSpline','EMG','CaSpike', ...
        'Gauss','doubleBoltzmann','Asym2Sigmoid'})} 
    options.t0 (1,1) {mustBeNumeric} = nan 
    options.bs (1,1) {mustBeNumeric} = nan
    options.x_ups (:,1) {mustBeNumeric} = linspace(x(1),x(end),numel(x)*10)
    options.mProf (:,1) logical = true(size(x))
    options.splineOrder (1,1) {mustBeInteger, mustBePositive} = 4
    options.spline_knots_seq (:,1) {mustBeNumeric, mustBeFinite} = 0
    options.addMoreWeightToEvent (1,1) logical = false
end

% start and end of event's mask
sE = find(options.mProf,1,'first');
eE = find(options.mProf,1,'last');
% x and y of event
x_E = x(sE:eE);
% y_E = y(sE:eE);
% x_E = x(sE:eE);
y_E = y;
y_E(~options.mProf) = -inf;
% max of event profile
[maxV, maxP] = max(y_E);
% step in x
dx = mean(diff(x));        
% options for lsqnonlin
optLsNonLin = optimoptions('lsqnonlin', ...
    'TolFun',1e-9, 'TolX',1e-9, 'TolCon',1e-9,...
    'MaxIter',1000, 'MaxFunEvals',3000, ...
    'Display','off');
% options for fitnlm
optFitNonLin = statset('TolFun',1e-9, 'TolX',1e-9,...
    'MaxIter',1000, 'MaxFunEvals',3000, ...
    'Display','off');
% first estimate t0, bs, tauR and tauD, used later to initialize fit
% fit rise and decay of event in time profile
% rise fun
funR = @(p,xR,yR) ((xR>=p(1)).*((1-exp(-(xR-p(1))./p(2))).*(p(3)) + p(4)) + ...
    (xR<p(1)).*(p(4)))-yR;
funR_eval = @(p,xR) (xR>=p(1)).*((1-exp(-(xR-p(1))./p(2))).*(p(3)) + p(4)) + ...
    (xR<p(1)).*(p(4));
% decay fun
funD = @(p,xD,yD) ( p(1).*exp(-(xD-p(2))/p(3)) + p(4) ) - yD;
funD_eval = @(p,xD) ( p(1).*exp(-(xD-p(2))/p(3)) + p(4) );
% rise part
xR = x(1:maxP);
yR = y(1:maxP);
% check for initial fit values
if isnan(options.t0) || ~isfinite(options.t0)
    % estimation of t0 and bs
    % baseline estimation
    t0FirstEstPx = sE - 1;
    if t0FirstEstPx <= 0
        t0FirstEstPx = 1; % first point
    end
    y0 = mean(yR(1:t0FirstEstPx));
    % amplitude est
    A = maxV - y0;
    try
        % estimate t0
        % find 25 and 75% from max of rise part
        p_75 = numel(yR) - ...
            find( (flipud(yR)-(y0 + (max(yR)-y0).*0.75))>0, 1, 'last') + 1;
        p_25 = numel(yR) - ...
            find( (flipud(yR)-(y0 + (max(yR)-y0).*0.25))<0, 1, 'first') + 1;
        if p_75 == p_25
            p_25 = p_75-1;
        end
        v_75 = yR(p_75);
        v_25 = yR(p_25);
        % construct line two-point form and get x at y=y0
        t0_est = xR(p_25)+( (y0-v_25)*(xR(p_75)-xR(p_25)) )/(v_75-v_25);
        tR_est = max((max(xR)-t0_est)*0.66, dx);
        if t0_est<xR(1), t0_est = xR(1); end
    catch
        t0_est = xR(maxP)-10;
        if t0_est<xR(1), t0_est = xR(1); end
        tR_est = dx;
    end
else
    % use the ones provided
    t0_est = options.t0;
    if t0_est < 0, t0_est = xR(1); end
    if isnan(options.bs) || ~isfinite(options.bs)
        y0 = mean(y(x<t0_est));
    else
        y0 = options.bs;
    end
    tR_est = max((max(xR)-t0_est)*0.66, dx);
    A = maxV - y0;
    % % correct for start of mask of event, also x start form 0
    % if options.t0 > maxP*dx - dx - 3
    %     p0(1) = maxP*dx - dx - 3;
    % end
end
% fit rise part of event
% p0 estimate [t0 tauR A bs]
p0Rise = [t0_est tR_est A y0];
try
    coeffRise = lsqnonlin(@(p)funR(p,xR,yR), p0Rise,...
        [min(xR) 0 0 min(yR)], ...
        [max(xR) max(xR)-min(xR) max(yR)-min(yR) max(yR)], ...
        optLsNonLin);
    t0 = coeffRise(1);
    bs = coeffRise(4);
    tauR_est = coeffRise(2);
    % check if tauR make sense
    if tauR_est > max(xR)-(t0-dx)
        tauR_est = (max(xR)-t0)*0.66;
    end
catch
    coeffRise = [];
    t0 = t0_est;
    bs = y0;
    tauR_est = tR_est;
end
%%%%%%%%%% estimate tauD
bs_end = mean(y(eE:end));
xD = x(maxP:end);
yD = y(maxP:end);
[~,p_tauD] = min( abs( yD - (0.33*(maxV-bs_end)+bs_end) ) );
tauD_est = p_tauD*dx;
% fit decay part
% [A, pos, tauD, bs]
p0Decay = [maxV-bs_end xD(1) tauD_est bs_end];
try
    coeffDecay = lsqnonlin(@(p)funD(p,xD,yD), p0Decay,...
        [0 min(xD) 0 min(yD)], ...
        [max(yD)-min(yD) max(xD) max(xD)-min(xD) max(yD)], ...
        optLsNonLin);
    tauD_est = coeffDecay(3);
catch
    coeffDecay = [];
end
% estimate amplitude
try
    A_est = maxV/(1-exp(-((maxP+sE-1)*dx-dx-t0)/tauR_est));
catch
    A_est = maxV-bs;
end

% setup model, p0 and bounds
switch model
    case '1expR1expD'
        % Lacampagne et al, 1999
        % coefficients: [t0, F01, tauR, A, t1, tauD, F02]
        p0 = [t0 bs tauR_est A_est x(maxP) tauD_est bs_end];
        lb = [min(x) min(y) 0 0 min(x) 0 min(y)];
        ub = [max(x) max(y) max(x)-min(x) max(y)-min(y) max(x) max(x)-min(x) max(y)];
        fitFun = @(p,x) (x<p(1)).*p(2) + ...
            (x>=p(1) & x<p(5)).*( p(2)+(1-exp(-(x-p(1))./p(3))).* (p(4)-p(2))) + ...
            (x>=p(5)).*( (p(4)-p(7)).*(1-exp(-(p(5)-p(1))./p(3))).*exp(-(x-p(5))./p(6))+p(7) );

    case 'spline'
        % setup number of knots for spline fitting
        if options.spline_knots_seq==0
            % every 5 points
            options.spline_knots_seq = linspace(x(1), x(end), ...
                2+fix(length(x)/5));
        end
        % coefficients: [spline coefficients]
        % splOrd = 4; % spline order, poly3
        spl0  = spap2(...
            augknt(options.spline_knots_seq, options.splineOrder), ...
            options.splineOrder, x, y);  % first guess of spline coefficients
        % newknt for a possibly better knot distribution
        try
            spl0  = spap2(newknt(spl0), options.splineOrder, x, y);
        catch
        end
        p0 = [spl0.coefs];
        fitFun = @(p,x) fnval(spmak(spl0.knots, spl0.coefs),x);

    case 'constThenSpline'
        if options.spline_knots_seq==0
            % every 10 ms
            options.spline_knots_seq = linspace(t0, x(end), ...
                2+fix(length(x(x>=t0))*dx/10));
        end
        % coefficients: [t0, F01, spline coefficients]
        % splOrd = 4; % spline order, poly3
        % first guess of spline coefficients
        spl0  = spap2(augknt(options.spline_knots_seq, options.splineOrder), ...
            options.splineOrder, x(x>t0), y(x>t0));
        % newknt for a possibly better knot distribution
        try
            spl0  = spap2(newknt(spl0), options.splineOrder, ...
                x(x>t0), y(x>t0));
        catch
        end
        p0 = [t0, bs, spl0.coefs];
        fitFun = @(p,x) (x <= p(1)) .*  p(2) + ...
            (x > p(1)) .* fnval(spmak(spl0.knots,p(3:end)),x);

    case 'EMG'
        % fit with exponentially modified gaussian
        % coefficients: [A m sd tau F0]
        sd_est = (sum(y>(bs+A_est/2))*dx)/2;
        A_est_EMG = (maxV-bs) / ...
            (exp(x(maxP)/tauD_est+sd_est^2/(2*tauD_est^2)-x(maxP)/tauD_est)* ...
            cdf('Normal',x(maxP),x(maxP)+(sd_est^2/tauD_est),sd_est));
        p0 = [A_est_EMG x(maxP) sd_est tauD_est bs];
        lb = [0 min(x) dx 0 min(y)];
        ub = [inf max(x) max(x)-min(x) max(x)-min(x) max(y)];
        fitFun = @(p,x) p(1).*exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - x./p(4)) .* ...
            cdf('Normal',x,p(2)+(p(3)^2/p(4)),p(3)) + p(5);

    case 'CaSpike'
        % Zahradnikova 2007, Kinetics of calcium spikes in rat cardiac myocytes
        % coefficients: [t0, F0, FM, tA, tI, FI]
        % coeff_n = {'t0','F01','Ampl','tauR','tauD','FI'};
        p0 = [t0 bs A_est tauR_est tauD_est 1];
        lb = [min(x) min(y) min(y) 0 0 0];
        ub = [max(x) max(y) max(y) max(x)-min(x) max(x)-min(x) 1];
        fitFun = @(p,x) (x < p(1)) .* p(2) + ...
            (x >= p(1)) .* (p(2) + (p(3)*((-3*p(4)*(2 + (-1 + p(6))*p(4))*(4*p(4) - 7*p(5))*(p(4) - 3*p(5)))./(2*exp((2*(x-p(1)))/p(4)))+...
            (3*p(4)*(p(4) - 3*p(5))*(p(4) - 2*p(5)))./exp((x-p(1))/p(4)) +(3*p(4)*(5 + (-1 + p(6))*p(4))*(2*p(4) - 5*p(5))*(p(4) - p(5)))./(5*exp((5*(x-p(1)))/p(4)))-...
            (p(4)*(6 + (-1 + p(6))*p(4))*(p(4) - 2*p(5))*(p(4) - p(5)))./(6*exp((6.*(x-p(1)))/p(4))) -(6*p(5)^3*(1 + (-1 + p(6))*p(5)))./exp((x-p(1))/p(5))+...
            (18*p(5)^3*(p(4) + p(5) + (-1 + p(6))*p(4)*p(5)))./(exp((p(4)^(-1) + p(5)^(-1)).*(x-p(1)))*(p(4) + p(5))) -...
            (18*p(5)^3*(p(4) + (2 + (-1 + p(6))*p(4))*p(5)))./(exp((2/p(4) + p(5)^(-1)).*(x-p(1)))*(p(4) + 2*p(5)))+...
            (6*p(5)^3*(p(4) + (3 + (-1 + p(6))*p(4))*p(5)))./(exp((3/p(4) + p(5)^(-1)).*(x-p(1)))* (p(4) + 3*p(5))) -...
            (3*p(4)*(4 + (-1 + p(6))*p(4))*(5*p(4)^2 - 20*p(4)*p(5) +17*p(5)^2))./(4.*exp((4.*(x-p(1)))/p(4))) +...
            (p(4)*(p(4)^2*(57 + 10*(-1 + p(6))*p(4)) - 3*p(4)*(84 + 13*(-1 + p(6))*p(4))*p(5) + (249 + 29*(-1 + p(6))*p(4))*p(5)^2))./(3.*exp((3.*(x-p(1)))/p(4)))-...
            ((-1 + p(6))*(p(4) - 3*p(5))*(p(4) - 2*p(5))*(p(4) - p(5))*(37*p(4)^4 + 252*p(4)^3*p(5) + 605*p(4)^2*p(5)^2 + 660*p(4)*p(5)^3 + 360*p(5)^4))./...
            (60*(p(4) + p(5))*(p(4) + 2*p(5))*(p(4) + 3*p(5))) +  (6*(-1 + p(6))*p(4)^2*(p(4) - 3*p(5))*(p(4) - 2*p(5)).*...
            cosh((x-p(1))/p(4)))./exp((2.*(x-p(1)))/p(4))))/...
            ((p(4) - 3*p(5))*(p(4) - 2*p(5))*(p(4) - p(5))));

    case 'Gauss'
        % fit with gaussian
        % coefficients: [F0, A, w, xc]
        % F0-baseline, A-amplitude, w-width of gaussian, xc-center
        A_est_gauss = (maxV-bs)*(sum(y>(bs+A_est/2))*dx)*sqrt(pi/(4*log(2)));
        p0 = [bs A_est_gauss sum(y>(bs+A_est/2))*dx x(maxP)];
        lb = [min(y) 0 dx min(x)];
        ub = [max(y) inf max(x)-min(x) max(x)];
        fitFun = @(p,x) p(1) + ...
            (p(2).*exp((-4.*log(2).*(x-p(4)).^2)./(p(3).^2))) ./ ...
            (p(3).*sqrt(pi/(4*log(2))));

    case 'doubleBoltzmann'
        % fit with double boltzmann function
        % coefficients: [F0, A, xH1, k1, xH2, k2]
        p0 = [bs A_est t0+tauR_est -tauR_est x(maxP)+tauD_est tauD_est];
        lb = [min(y) 0 min(x) -(max(x)-min(x)) min(x) 0 ];
        ub = [max(y) inf max(x) 0 max(x) max(x)-min(x)];
        fitFun = @(p,x) p(1) + p(2) .* ...
            ( 1./(1+exp((x-p(3))./p(4))) ) .* ...
            ( 1./(1+exp((x-p(5))./p(6))) ) ;

    case 'Asym2Sigmoid'
        % originLab
        % coefficients:[F0, A, xc, w1, w2, w3]
        % Meanings: y0 = offset, xc = center, A = amplitude,
        % w1 = full width of half maximum,
        % w2 = variance of low-energy side,
        % w3 = variance of high-energy side.
        % Lower Bounds: w1 > 0.0, w2 > 0.0, w3 > 0.0
        p0 = [bs A_est x(maxP) sum(y>(bs+A_est/2))*dx tauR_est tauD_est];
        lb = [min(y) 0 min(x) 0 0 0];
        ub = [max(y) inf max(x) max(x)-min(x) max(x)-min(x) max(x)-min(x)];
        fitFun = @(p,x) p(1) + p(2) .* ...
            ( 1./(1+exp(-(x-p(3)+p(4)/2)./p(5))) ) .* ...
            ( 1 - 1./(1+exp(-(x-p(3)-p(4)/2)./p(6))) ) ;
end

% observation weights
wE = ones(size(y));
if options.addMoreWeightToEvent
    % put more weights to event than to baseline
    wE(options.mProf) = (1 + y(options.mProf)./max(y(options.mProf))).^2;
end
fitFunDiff = @(p,x,y) wE.*(y-fitFun(p,x));
% fit whole event
switch  model
    case 'spline'
        coef = p0;
        tauD = tauD_est;
        tauR = tauR_est;
    case 'constThenSpline'
        coef = p0;
        tauD = tauD_est;
        tauR = tauR_est;
        bs = coef(2);
        t0 = coef(1);
    otherwise
        try
            % fit profile with model
            try
                % no parameters constraints
                mdl = fitnlm(x, y, fitFun, p0, ...
                    'Weights',wE, 'Options',optFitNonLin);
                % check if fit is better than constant 
                assert(mdl.ModelFitVsNullModel.Pvalue<0.05)
                coef = mdl.Coefficients.Estimate;
            catch
                coef = lsqnonlin(@(p)fitFunDiff(p,x,y), ...
                    p0, lb, ub, optLsNonLin);
            end
            % check if there are no Nan or inf once function is evaluated
            assert(all(isfinite(fitFun(coef,x))))
            % assign t0 and baseline parameters form fits with
            % different models
            switch model
                case '1expR1expD'
                    % [t0, F01, tauR, A, t1, tauD, F02]
                    bs = coef(2);
                    t0 = coef(1);
                    tauR = coef(3);
                    tauD = coef(6);
                case 'EMG'
                    % [A m sd tau F0]
                    bs = coef(5);
                    tauR = tauR_est;
                    tauD = coef(4);
                case 'CaSpike'
                    % [t0, F0, FM, tA, tI, FI]
                    bs = coef(2);
                    t0 = coef(1);
                    tauR = coef(4);
                    tauD = coef(5);
                case 'Gauss'
                    % [F0, A, w, xc]
                    bs = coef(1);
                    tauR = tauR_est;
                    tauD = tauD_est;
                case 'doubleBoltzmann'
                    % [F0, A, xH1, k1, xH2, k2]
                    bs = coef(1);
                    tauR = tauR_est;
                    tauD = tauD_est;
                case 'Asym2Sigmoid'
                    % [F0, A, xc, w1, w2, w3]
                    bs = coef(1);
                    tauR = tauR_est;
                    tauD = tauD_est;
            end
        catch
            % if fails fit with spline
            model = 'spline';
            % setup number of knots for spline fitting
            if options.spline_knots_seq==0
                % every 5 points
                options.spline_knots_seq = linspace(x(1), x(end), ...
                    2+fix(length(x)/5));
            end
            % coefficients: [spline coefficients]
            % splOrd = 4; % spline order, poly3
            spl0  = spap2(...
                augknt(options.spline_knots_seq, options.splineOrder), ...
                options.splineOrder, x, y);  % first guess of spline coefficients
            % newknt for a possibly better knot distribution
            try
                spl0  = spap2(newknt(spl0), options.splineOrder, x, y);
            catch
            end
            p0 = [spl0.coefs];
            fitFun = @(p,x) fnval(spmak(spl0.knots, spl0.coefs),x);
            coef = p0;
            tauD = tauD_est;
            tauR = tauR_est;
        end
end

% get fit
yFit = fitFun(coef, x);
try
    switch model
        case {'1expR1expD', 'constThenSpline'}
            yFit_ups = interp1(x, yFit, options.x_ups, 'makima');
        otherwise
            yFit_ups = fitFun(coef,options.x_ups);
    end
catch
    options.x_ups = [];
    yFit_ups = [];
end

% figure('Name',sprintf("%s",model))
% plot(x,y,'ok')
% hold on
% plot(options.x_ups,yFit_ups,'or')
% plot(options.x_ups,yFit_ups2,'ob')
% 
% plot(x,fitFun(p0,x),'g')
% plot(x,fitFun(coef,x),'b')
                  
% create output structure
out = struct('fitModel',model,...
    'coeff',coef,...
    'bs',bs,...
    't0',t0,...
    'tauD',tauD,...
    'tauR',tauR,...
    'x',x, ...
    'y',y,...
    'evntMaxV',maxV, ...
    'evntMaxPosPx',maxP, ...
    'yFit',yFit(:),...
    'x_ups',options.x_ups(:),...
    'yFit_ups',yFit_ups(:), ...
    'yFitFunR',funR_eval, ...
    'coeffRise',coeffRise,...
    'yFitFunD',funD_eval, ...
    'coeffDecay',coeffDecay);
end
