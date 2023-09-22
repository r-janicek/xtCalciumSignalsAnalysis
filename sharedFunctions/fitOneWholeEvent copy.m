function out = fitOneWholeEvent(p0,x,x_ups,y,piecewise,model,param1,param2)

% parameters of function
% p0 = coeficients
% x = independent variable
% x_ups = upscaled
% y = data to fit
% piecewise = options yes or no
% models:
% 1. 1expR1expD (single exponential rise and decay)
% 2. spline
% 3. CaSpikeFun
% 4. EMG (exponentially modified gaussian)
% 5. Asym2Sigmoid
% 6. doubleBoltzmann

% write function where it will store all fitting funcions, output will be
% function, number of coeff, coeff names. This then will be used for 
% fitting single and multiple peaks fitting function within whole analysis
% also incorporate sign of fit and also baseline
keyboard
% data
x = x(:);
y = y(:);
if ~isempty(param1)
    nK = param1;
    splOrd = param2;
end
% step in x
dx = mean(diff(x));
        
% options for fmincon
opt = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
    'MaxIter',1000,'MaxFunEvals',3000);

% do fitting
switch piecewise
    
    case 'yes'
        
        % fit rise of spark in time profile
        optFit = optimoptions('lsqnonlin','TolFun',1e-9,'TolX',1e-9,'MaxIter',1000,'MaxFunEvals',3000,'Display','off');
       
        funR = @(p,xR,yR) ((xR>=p(1)).*((1-exp(-(xR-p(1))./p(2))).*(p(3)) + p(4)) + ...
            (xR<p(1)).*(p(4)))-yR;
        
        funR_eval = @(p,xR) (xR>=p(1)).*((1-exp(-(xR-p(1))./p(2))).*(p(3)) + p(4)) + ...
             (xR<p(1)).*(p(4));
               
        funD = @(p,xD,yD) ( p(1).*exp(-(xD-p(2))/p(3)) + p(4) ) - yD; 
        funD_eval = @(p,xD) ( p(1).*exp(-(xD-p(2))/p(3)) + p(4) ); 
                
        [maxV, maxP] = max(y);
         
        xR = x(1:maxP);
        yR = y(1:maxP);
        
        if p0(1) > maxP*dx - 3           
            p0(1) = maxP*dx - 5;
        end
        
        if p0(1) < 0, p0(1) = 1; end
        
        p0Rise = [p0(1) 5 maxV p0(2)];  % [t0 tauR A bs]
        % fit rise 
     
        try 
            coeffRise = lsqnonlin(@(p)funR(p,xR,yR),p0Rise,[zeros(1,3),min(yR)],ones(1,4).*inf,optFit);
            t0 = coeffRise(1);
            bs = coeffRise(4);
        
        catch 
            t0 = p0(1);
            bs = p0(2);          
        end  
        
        bs_end = mean(y(end-2:end)); 
        
        % estimate tauR     
        t0_p = floor(t0/dx); % t0 in points
        if t0_p<1, t0_p = 1; end
        
        [~,p_tauR] = min( abs( yR(t0_p:end) - (0.66*(maxV-bs)+bs) ) );                
        tauR_est = p_tauR*dx;
              
        % estimate tauD
        xD = x(maxP+1:end);
        yD = y(maxP+1:end);
        [~,p_tauD] = min( abs( yD - (0.33*(maxV-bs_end)+bs_end) ) );   
        tauD_est = p_tauD*dx;
        
        % fit decay
        try
            p0Decay = [maxV xD(1) tauD_est bs_end];
            coeffDecay = lsqnonlin(@(p)funD(p,xD,yD),p0Decay,[zeros(1,3),min(yR)],ones(1,4).*inf,optFit);
            tauDFitExp = coeffDecay(3);
        catch
            tauDFitExp = tauD_est;
        end
        
        % estimate ampl 
        A_est = maxV/(1-exp(-(maxP*dx-t0)/tauR_est));
        %keyboard                    
%          figure
%          plot(x,y,'ob')
%          hold on         
%          plot(xR,yR,'ok')
%          plot(xD,yD,'or')
%          plot(xR,funR_eval(coeffRise,xR),'k','LineWidth',2)
%          plot(xD,funD_eval(coeffDecay,xD),'r','LineWidth',2)
                  
        p0(1) = t0;
        p0(2) = bs;
       % tauD = 5
                      
        switch model
                       
            case '1expR1expD'
                % Lacampagne et al, 1999
                % coefficients: [t0, F01, tauR, A, t1, tauD, F02]
                fitFun = @(p,x) (x<p(1)).*p(2) + ...
                                (x>=p(1) & x<=p(5)).*( p(2)+(1-exp(-(x-p(1))./p(3))).*p(4) ) + ...
                                (x>p(5)).*( p(4).*(1-exp(-(p(5)-p(1))./p(3))).*exp(-(x-p(5))./p(6))+p(7) );     
                p0(3) = tauR_est;
                p0(4) = A_est;
                p0(6) = tauD_est;
                
                coeff_n = {'t0','F01','tauR','Ampl','tmax','tauD','F02'};
                                                       
             
            case '1expR1expD_2'
                % Lacampagne et al, 1999
                % coefficients: [t0, F01, tauR, A, t1, tauD, F02]
                 
                p0(3) = tauR_est;
                p0(4) = A_est;
                p0(6) = tauD_est;
                
                coeff_n = {'t0','F01','tauR','Ampl','tmax','tauD','F02'};
                           
                fitFun = @(p,x) (x<p(1)).*p(2) + ...
                                (x>=p(1) & x<p(5)).*( p(2)+(1-exp(-(x-p(1))./p(3))).* (p(4)-p(2))) + ...
                                (x>=p(5)).*( (p(4)-p(7)).*(1-exp(-(p(5)-p(1))./p(3))).*exp(-(x-p(5))./p(6))+p(7) );  
                            
                     
            case '1expR1expD_3'
                % Lacampagne et al, 1999
                % coefficients: [t0, F01, tauR, A, t1, tauD, F02]
                 
                p0(3) = tauR_est;
                p0(4) = A_est;
                p0(6) = tauD_est;
                
                coeff_n = {'t0','F01','tauR','Ampl','tmax','tauD','F02'};
                                                      
                fitFun = @(p,x) (x<p(1)).* p(2) + ...
                                (x>=p(1)).* ( p(2) + (1-exp(-(x-p(1))./p(3))) .* p(4) .* (exp(-(x-p(5))./p(6))+p(7) ) );        
                            
                            
            case 'spline'
                % coefficients: [t0, F01, spline coefficients]
                %nK = ceil((max( x(x>p0(1)))-min( x(x>p0(1))))/5);
                %splOrd = 4; % spline order, poly3
                spl0  = spap2(nK, splOrd, x(x>p0(1)), y(x>p0(1)));  % first guess of spline coefficients
                % newknt for a possibly better knot distribution
                knots = newknt(spl0);
                spl0  = spap2(knots, splOrd, x(x>p0(1)), y(x>p0(1)));
             
                nK = numel(spl0.knots);

                p0 = [p0(1), p0(2), spl0.coefs];
                coeff_n = [{'t0','F01'}, repmat({'splCoef'},size(spl0.coefs))];
                
                fitFun = @(p,x) (x <= p(1)) .*  p(2) + (x > p(1)) .* fnval(spmak(spl0.knots,p(3:end)),x);
               

            case 'CaSpikeFun'
                %coefficients: [t0, F0, FM, tA, tI, FI]
                p0(3) = A_est;
                p0(4) = tauR_est;
                p0(5) = tauD_est;

                coeff_n = {'t0','F01','Ampl','tauR','tauD','FI'};
                
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
                
                
        end
        
        % put more weights to event than to baseline
        wE = ones(size(y));
        wE(x>p0(1)) = ( 1 + y(x>p0(1))./max(y(x>p0(1))) ).^2;
        
        % do fitting
        % create function which will be optimize, return scalar
        fitFunSum = @(p,x,y) sum( wE .* (y-fitFun(p,x)).^2 );
        
        switch  model
            
            case 'spline'

                coef = p0;
                tauD = tauDFitExp;   
               
                
            case {'1expR1expD','CaSpikeFun'}
              
                try
                    lb = zeros(size(p0));
                    ub = inf(size(p0));
                    
%                      figure
%                      plot(x,y)
%                      hold on
%                      plot(x,fitFun(p0,x),'r')
%                     plot(x,fitFun(coef,x),'g')
                    
                    switch  model
                        case '1expR1expD'
                            % first fit with fixed t0, t1, baseline and tauR = these were estimated
                            lb(1) = p0(1); ub(1) = p0(1);
                            lb(2) = p0(2); ub(2) = p0(2);
                            lb(3) = p0(3); ub(3) = p0(3);
                            lb(5) = p0(5); ub(5) = p0(5);
                        case 'CaSpikeFun'
                            lb = [];
                            ub = [];
                    end
                    coef = fmincon(@(p)fitFunSum(p,x,y),p0,[],[],[],[],lb,ub,[],opt);
                    
                    % refit with no bound, no weights
                    lb = zeros(size(p0));
                    ub = inf(size(p0));
                    wE = ones(size(y));
                    fitFunSum = @(p,x,y) sum( wE .* (y-fitFun(p,x)).^2 );
                    coef = fmincon(@(p)fitFunSum(p,x,y),coef,[],[],[],[],lb,[],[],opt);
                    
                    switch  model
                        case '1expR1expD'
                            tauD = coef(6);
                            
                        case 'CaSpikeFun'
                            tauD = coef(5);
                    end
                    
                catch
                    % refit with spline
               
                    % coefficients: [t0, F01, spline coefficients]
                    nK = ceil((max( x(x>p0(1)))-min( x(x>p0(1))))/5);
                    splOrd = 4; % spline order, poly3
                    spl0  = spap2(nK, splOrd, x(x>p0(1)), y(x>p0(1)));  % first guess of spline coefficients
                    % newknt for a possibly better knot distribution
                    knots = newknt(spl0);
                    spl0  = spap2(knots, splOrd, x(x>p0(1)), y(x>p0(1)));
                    
                    nK = numel(spl0.knots);
                    
                    p0 = [p0(1), p0(2), spl0.coefs];
                    
                    fitFun = @(p,x) (x <= p(1)) .*  p(2) + (x > p(1)) .* fnval(spmak(spl0.knots,p(3:end)),x);
                    
                    coef = p0;
                    coeff_n = [{'t0','F01'}, repmat({'splCoef'},size(spl0.coefs))];
                    model = 'spline';
                    
                    tauD = tauDFitExp;
                    
                end
                
        end
        
        yFit = fitFun(coef,x);
        if ~isempty(x_ups)
            yFit_ups = pchip(x,yFit,x_ups);
        else
            x_ups = [];
            yFit_ups = [];
        end
        
        bs = coef(2);
        t0 = coef(1);
        
                    
    case 'no'
        
        switch  model
            
            case 'expModGauss'               
                % fit with exponentially modified gaussian
                % coefficients: [A,m,sd,tau,F0]
                fitFun = @(p,x) p(1).*exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - x./p(4)).*cdf('Normal',x,p(2)+(p(3)^2/p(4)),p(3)) + p(5);
                
                fitFunSum = @(p,x,y) sum( (y-fitFun(p,x)).^2 );
                coef = fmincon(@(p)fitFunSum(p,x,y),p0,[],[],[],[],zeros(size(p0)),[],[],opt);
                
                bs = coef(5);
                tauD = coef(4);
                
                coeff_n = {'A','m','sd','tauD', 'F01'};
             
                
            case 'doubleBoltzmann'
                % fit with double boltzmann function
                % coefficients: [F0, A, xH1, k1, xH2, k2]
                fitFun = @(p,x) p(1) + p(2) .* ( 1./(1+exp((x-p(3))./p(4))) ) .* ( 1./(1+exp((x-p(5))./p(6))) ) ;
                
                fitFunSum = @(p,x,y) sum( (y-fitFun(p,x)).^2 );
                coef = fmincon(@(p)fitFunSum(p,x,y),p0,[],[],[],[],[0 0 0 0 0 -inf],[],[],opt);
                
                bs = coef(1);
                tauD = nan;
                
                coeff_n = {'F0','A','xH1','k1', 'xH2','k2'};
                         
                    
            case 'Asym2Sigmoid'
                % originLab
                % coefficients:[F0, A, xc, w1, w2, w3]
                % Meanings: y0 = offset, xc = center, A = amplitude, w1 = full width of half maximum,
                % w2 = variance of low-energy side, w3 = variance of high-energy side.
                % Lower Bounds: w1 > 0.0, w2 > 0.0, w3 > 0.0
                fitFun = @(p,x) p(1) + p(2) .* ( 1./(1+exp(-(x-p(3)+p(4)/2)./p(5))) ) .* ( 1 - 1./(1+exp(-(x-p(3)-p(4)/2)./p(6))) ) ;

                fitFunSum = @(p,x,y) sum( (y-fitFun(p,x)).^2 );
                coef = fmincon(@(p)fitFunSum(p,x,y),p0,[],[],[],[],zeros(size(p0)),[],[],opt);
                
                bs = coef(1);
                tauD = nan;
                
                coeff_n = {'F0','A','xc','w1','w2', 'w3'};
            
                
            case 'spline'
                % fit with spline
               
                %nK = ceil( (max(x)-min(x))/1 ); % fit every 1 um
                %splOrd = 3; % spline order, poly3
                spl0  = spap2(nK, splOrd, x, y);  % first guess of spline coefficients
                % newknt for a possibly better knot distribution
                knots = newknt(spl0);
                spl0  = spap2(knots, splOrd, x, y);
                
                fitFun = @(p,x) fnval(spmak(spl0.knots,spl0.coefs),x);
                coef = spl0.coefs;
                coeff_n = repmat({'splCoef'},size(spl0.coefs));
                bs = nan;
                tauD = nan;
                
        end
        
        yFit = fitFun(coef,x);
        if ~isempty(x_ups)
            yFit_ups = fitFun(coef,x_ups);
        else
            x_ups = [];
            yFit_ups = [];
        end
       
        t0 = nan;
        
end

% create output structure
% out = struct('fitModel',model,...
%     'coeff',coef,...
%     'bs',bs,...
%     't0',t0,...
%     'tauD',tauD,...
%     'x',x,'y',y,...
%     'yFit',yFit(:),...
%     'x_ups',x_ups(:),...
%     'yFit_ups',yFit_ups(:));

% same structure as from multiple peak fitting
out.t = x;
out.wholeFit = yFit(:);
out.allEventsFit = yFit(:);
out.baselineFit = zeros(size(yFit(:)))+bs;
out.individualEventsFits = [{'event 1'}; num2cell(yFit(:))];
out.coefficientsOfFittedEvents = [[{''},coeff_n];[{'event 1'},num2cell(coef)]];
out.eventModel = model;
out.bs = bs;
out.t0 = t0;

if exist('x_ups', 'var') && ~isempty(x_ups)
    out.t_ups.wholeFit = yFit_ups(:);
    out.t_ups.allEventsFit = yFit_ups(:);
    out.t_ups.baselineFit = zeros(size(yFit_ups(:)))+bs;
    out.t_ups.individualEventsFits = [ {'event 1'} ; num2cell(yFit_ups(:)) ];
    out.t_ups.t_ups = x_ups(:);
end



