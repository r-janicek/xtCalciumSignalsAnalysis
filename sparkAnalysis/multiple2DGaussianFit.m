function out = multiple2DGaussianFit(...
    T,X,T_ups,X_ups,imgPair,imgPair_m,imgPair_events_m,...
    prevFitParams,sPP_pair,...
    fitFun,bsWholeProfFit)

% T,X,D = matrices , time, spatial and image
% _ups = upscaled
% maxOfEvents = peaks value
% posOfEvents_t = position of peak in t
% posOfEvents_x = position of peak in x
% sd_t,expTauD_t,sPP_pair,bsWholeProfFit

N_events = numel(prevFitParams.maxOfEvents);
pxSzT = mean(diff(T(1,:)));
pxSzX = mean(diff(X(:,1)));

% event fit function
[fun_event,p0,A_ineq,b_ineq,lb,ub,n_coef_event,coeff_n] = ...
    setFitFunction(fitFun,N_events,prevFitParams, ...
    sPP_pair,bsWholeProfFit,pxSzT,pxSzX,...
    T,X,imgPair,imgPair_events_m);

% do fitting
events_coeff = ...
    fitOfData(p0,T,X,imgPair,imgPair_m,n_coef_event,fun_event,N_events,...
    A_ineq,b_ineq,lb,ub);

% calculate final fit and fits of events 
[F,F_events,F_individualEvents] = ...
    sumEvents(events_coeff,T,X,n_coef_event,fun_event,N_events);

% upscaled
[F_ups,F_events_ups,F_individualEvents_ups] = ...
    sumEvents(events_coeff,T_ups,X_ups,n_coef_event,fun_event,N_events);

% create final output structure
events_n = num2cell(zeros(1,N_events));
coefficientsOfEvents = zeros(N_events,n_coef_event);
individualEventsFits = cell([2,numel(events_n)]);
individualEventsFitsUps = cell([2,numel(events_n)]);
for i=1:N_events
    events_n(1,i) = {sprintf('event %d',i)};
    coefficientsOfEvents(i,:) = events_coeff( 1+(i-1)*n_coef_event : n_coef_event+(i-1)*n_coef_event );
    individualEventsFits(:,i) = [ events_n(1,i) ; {F_individualEvents(:,:,i)} ];
    individualEventsFitsUps(:,i) = [ events_n(1,i) ; {F_individualEvents_ups(:,:,i)} ];
end

out.T = T;
out.X = X;
out.img = imgPair;
out.wholeFit = F;
out.allEventsFit = F_events;
out.baselineEst = events_coeff(end);
out.individualEventsFits = individualEventsFits;
out.coefficientsOfFittedEvents = [[{''},coeff_n];[events_n',num2cell(coefficientsOfEvents)]];
out.eventModelName = fitFun;
out.eventModelFun = fun_event;

% upscaled
out.ups.T_ups = T_ups;
out.ups.X_ups = X_ups;
out.ups.wholeFit = F_ups;
out.ups.allEventsFit = F_events_ups;
out.ups.individualEventsFits = individualEventsFitsUps;

end


%%%%%%%%%%%%%%%%%%%%%%%% functions for fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do fitting
function events_coeff = fitOfData( ... 
    p0,T,X,imgPair,imgPair_m,...
    n_coef_event,fun_event,N_events,A_ineq,b_ineq,lb,ub)


% Linear inequality constraints = A_ineq,b_ineq
% lower and upper bounds = lb, ub

% global
% tic
% %%%%%%%%% do fitting %%%%%%%%%%%
% %fit with weights to push to fit events data preferentially
% W = ones(size(imgPair));
% % W(imgPair_m) = ( 1 + imgPair(imgPair_m)./max(imgPair(imgPair_m)) ).^1;
% % global optimalization
% options = optimoptions(@simulannealbnd,'MaxIterations',400,...
%     'HybridFcn','fmincon','MaxTime',300);
% events_coeff = simulannealbnd(...
%     @(p)fitFunSumOfSquaersW(p,T,X,imgPair,n_coef_event,fun_event,N_events,W),...
%     p0,lb,ub,options);
% toc
% 
% 
% events_coeff = patternsearch(...
%     @(p)fitFunSumOfSquaersW(p,T,X,imgPair,n_coef_event,fun_event,N_events,W),...
%     p0,[],[],[],[],lb,ub);
%  toc

tic
% local 
% FMINCON
options = optimoptions('fmincon','TolFun',1e-12,'TolX',1e-12,...
    'TolCon',1e-12,...
    'MaxIter',1000,'MaxFunEvals',3000,...
    'UseParallel',false,'Display','none',...
    'OutputFcn','');

% % fit with weights to push to fit events data preferentially, lock baseline
% W = zeros(size(imgPair));
% W(imgPair_m) = ( 1 + imgPair(imgPair_m)./max(imgPair(imgPair_m)) ).^4;
%  
% % fit first
% lb_singleEvnt = lb;
% lb_singleEvnt(9:end-1) = p0(9:end-1);
% ub_singleEvnt = ub;
% ub_singleEvnt(9:end-1) = p0(9:end-1);
%  
% events_coeff = ...
%     fmincon(@(p)fitFunSumOfSquaersW(p,T,X,imgPair,n_coef_event,fun_event,N_events,W),...
%     p0,A_ineq,b_ineq,[],[],lb_singleEvnt,ub_singleEvnt,[],options);
% 
% % fit second
% lb_singleEvnt = lb;
% lb_singleEvnt(1:8) = events_coeff(1:8);
% lb_singleEvnt(end) = events_coeff(end);
% ub_singleEvnt = ub;
% ub_singleEvnt(1:8) = events_coeff(1:8);
% ub_singleEvnt(end) = events_coeff(end);
% 
% events_coeff = ...
%     fmincon(@(p)fitFunSumOfSquaersW(p,T,X,imgPair,n_coef_event,fun_event,N_events,W),...
%     events_coeff,A_ineq,b_ineq,[],[],lb_singleEvnt,ub_singleEvnt,[],options);

% refit, to get better baseline estimation, lock all other parameters
% except baseline, take only beggining of image, where baseline sjould be
% real
W = ones(size(imgPair));
% start of pair of events
ind_E_s = find(sum(imgPair_m,1)>0,1,'first');
W(:,ind_E_s:end) = 0;

lb_bs = [p0(1:end-1),min(imgPair(:))];
ub_bs = [p0(1:end-1),max(imgPair(:))];
events_coeff = ...
    fmincon(@(p)fitFunSumOfSquaersW(p,T,X,imgPair,n_coef_event,fun_event,N_events,W),...
    p0,A_ineq,b_ineq,[],[],lb_bs,ub_bs,[],options);

% fit both events together
% refit with weights, lock baseline, 
W = zeros(size(imgPair));
W(imgPair_m) = ( 1 + imgPair(imgPair_m)./max(imgPair(imgPair_m)) ).^4;
lb(end) = events_coeff(end);
ub(end) = events_coeff(end);
events_coeff = ...
    fmincon(@(p)fitFunSumOfSquaersW(p,T,X,imgPair,n_coef_event,fun_event,N_events,W),...
    events_coeff,A_ineq,b_ineq,[],[],lb,ub,[],options);

toc 
 



% %events_coeff = p0;
% 
% [F,F_events,F_individualEvents] = ...
%     sumEvents(events_coeff,T,X,n_coef_event,fun_event,N_events);
% 
% figure
% mesh(T,X,imgPair,'LineStyle','-',...
%     'LineWidth',0.5,'FaceColor','none',...
%     'EdgeColor','k','EdgeAlpha',0.4)
% hold on
% mesh(T,X,F_individualEvents(:,:,1)+events_coeff(end),'LineStyle','-',...
%     'LineWidth',0.5,'FaceColor','b',...
%     'EdgeColor','none','EdgeAlpha',0.1,'FaceAlpha',0.4)
% mesh(T,X,F_individualEvents(:,:,2)+events_coeff(end),'LineStyle','-',...
%     'LineWidth',0.5,'FaceColor','g',...
%     'EdgeColor','none','EdgeAlpha',0.1,'FaceAlpha',0.4)
% mesh(T,X,F,'LineStyle','-',...
%     'LineWidth',0.5,'FaceColor','k',...
%     'EdgeColor','none','EdgeAlpha',0.1,'FaceAlpha',0.2)
% % % 
% %surf(T,X,F,'FaceAlpha',0.6,'EdgeColor','none','FaceColor','interp')

       

end


% set fit function
function [fun_event,p0,A_ineq,b_ineq,lb,ub,n_coef_event,coeff_n] = ...
    setFitFunction(fitFun,N_events,...
    prevFitParams,sPP_pair,bs,pxSzT,pxSzX,...
    T,X,imgPair,imgPair_events_m)

% estimates of bounds for parameters
sE = zeros(1,N_events);
eE = zeros(size(sE));
lb_m1 = zeros(size(sE));
ub_m1 = zeros(size(sE));
lb_m2 = zeros(size(sE));
ub_m2 = zeros(size(sE));
for i=1:N_events
    % start and end of event in x direction (time)
    [~,peakPosX] = min(abs(X(:,1)-prevFitParams.posOfEvents_x(i)));
    sE(i) = find(imgPair_events_m(peakPosX,:,i),1,'first');
    eE(i) = find(imgPair_events_m(peakPosX,:,i),1,'last');
    
    % bounds for positions parameters
    lb_m1(i) = sE(i)*pxSzT-pxSzT; if lb_m1(i)<0, lb_m1(i)=0; end
    ub_m1(i) = eE(i)*pxSzT;
    lb_m2(i) = min(min(X(imgPair_events_m(:,:,i))));
    ub_m2(i) = max(max(X(imgPair_events_m(:,:,i))));  
end


switch fitFun
    
    case '2D_EMGt'
        % 2D exponentially modified gaussian in time direction only
        % coefficients: [A, peakPos_t, sd_t, expTauD_t, peakPos_x, sd_x]       
        
        fun_event = @(p,x,y) p(1) .* exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - x./p(4)) .* ...
            cdf('Normal',x,p(2)+(p(3)^2/p(4)),p(3)) .* ...
            exp(-((y-p(5)).^2./(2*p(6).^2)));
        
        n_coef_event = 6;
        coeff_n = {'A', 'peakPos_t', 'sd_t', 'expTauD_t', 'peakPos_x', 'sd_x', 'F0'};

        p0 = zeros(1,N_events*n_coef_event);
        lb = zeros(1,N_events*n_coef_event);
        ub = inf(1,N_events*n_coef_event);
          A_ineq = zeros(4*N_events,N_events*n_coef_event+1);
        b_ineq = zeros(size(A_ineq,1),1);
keyboard
        % coefficients: [A, peakPos_t, sd_t, expTauD_t, peakPos_x, sd_x]
        for i=1:N_events
            
            A = 2*pi*prevFitParams.sd_t(i)*...
                prevFitParams.sd_x(i)*...
                prevFitParams.maxOfEvents(i);
            p0(1 + (i-1)*n_coef_event) = A;
            p0(2 + (i-1)*n_coef_event) = prevFitParams.posOfEvents_t(i);
            p0(3 + (i-1)*n_coef_event) = prevFitParams.sd_t(i);
            p0(4 + (i-1)*n_coef_event) = prevFitParams.expTauD_t(i);
            p0(5 + (i-1)*n_coef_event) = prevFitParams.posOfEvents_x(i);
            p0(6 + (i-1)*n_coef_event) = prevFitParams.sd_x(i);
            
            lb(2 + (i-1)*n_coef_event) = min(sPP_pair(i),lb_m1(i));
            try
                ub(2 + (i-1)*n_coef_event) = min(sPP_pair(i+1)-pxSzT,ub_m1(i));
            catch
                ub(2 + (i-1)*n_coef_event) = ub_m1(i);
            end
            
            lb(3 + (i-1)*n_coef_event) = pxSzT;
            ub(3 + (i-1)*n_coef_event) = ub_m1(i)-lb_m1(i);
            lb(4 + (i-1)*n_coef_event) = pxSzT;
            ub(4 + (i-1)*n_coef_event) = ub_m1(i)-lb_m1(i);
            lb(5 + (i-1)*n_coef_event) = lb_m2(i);
            ub(5 + (i-1)*n_coef_event) = ub_m2(i);
            lb(6 + (i-1)*n_coef_event) = pxSzX;
            ub(6 + (i-1)*n_coef_event) = ub_m2(i)-lb_m2(i);
            
             % Linear inequality constraints
            % sample mean of expModGaussian should be in image (m=mu+tau)
            A_ineq(2*i-1,6 + (i-1)*n_coef_event) = 1;
            A_ineq(2*i-1,8 + (i-1)*n_coef_event) = 1;
            A_ineq(2*i,6 + (i-1)*n_coef_event) = -1;
            A_ineq(2*i,8 + (i-1)*n_coef_event) = -1;
            b_ineq(2*i-1) = ub(6 + (i-1)*n_coef_event);
            b_ineq(2*i) = -lb(6 + (i-1)*n_coef_event);
            
        end
        % baseline
        lb(end) = min(imgPair(:));
        ub(end) = max(imgPair(:));
        
        % add baseline estimate at the end
        p0 = [p0,bs];
        
                 
    case '2D_EMGxt'

        % 2D exponentially modified gaussian in each direction
        % coefficients: [A, peakPos_t, sd_t, expTauD_t, peakPos_x, sd_x, expTauD_x, F0]
        % p0 = [2*pi*sd_t*sd_y*A m1 sd_t expTauD_t m2 sd_x expTauD_x];
        
        fun_event = @(p,x,y) p(1) .* ...
            exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - x./p(4)) .* ...
            cdf('Normal',x,p(2)+(p(3)^2/p(4)),p(3)) .* ...     
            exp(p(5)/p(7) + p(6)^2/(2*p(7)^2) - y./p(7)) .* ...
            cdf('Normal',y,p(5)+(p(6)^2/p(7)),p(6));
        
        n_coef_event = 7;
        coeff_n = {'A', 'peakPos_t', 'sd_t', 'expTauD_t', ...
            'peakPos_x', 'sd_x', 'expTauD_x'};

        p0 = zeros(1,N_events*n_coef_event);
        lb = zeros(1,N_events*n_coef_event+1);
        ub = inf(1,N_events*n_coef_event+1);
        A_ineq = zeros(4*N_events,N_events*n_coef_event+1);
        b_ineq = zeros(size(A_ineq,1),1);
keyboard
        for i=1:N_events
            
            A = 2*pi*prevFitParams.sd_t(i)*...
                prevFitParams.sd_x(i)*...
                prevFitParams.maxOfEvents(i);
        
            p0(1 + (i-1)*n_coef_event) = A;
            p0(2 + (i-1)*n_coef_event) = prevFitParams.posOfEvents_t(i);
            p0(3 + (i-1)*n_coef_event) = prevFitParams.sd_t(i);
            p0(4 + (i-1)*n_coef_event) = prevFitParams.expTauD_t(i);
            p0(5 + (i-1)*n_coef_event) = prevFitParams.posOfEvents_x(i);
            p0(6 + (i-1)*n_coef_event) = prevFitParams.sd_x(i);
            p0(7 + (i-1)*n_coef_event) = prevFitParams.expTauD_x(i);
            
            lb(2 + (i-1)*n_coef_event) = min(sPP_pair(i),lb_m1(i));
            try
                ub(2 + (i-1)*n_coef_event) = min(sPP_pair(i+1)-pxSzT,ub_m1(i));
            catch
                ub(2 + (i-1)*n_coef_event) = ub_m1(i);
            end
            
            lb(3 + (i-1)*n_coef_event) = pxSzT;
            ub(3 + (i-1)*n_coef_event) = ub_m1(i)-lb_m1(i);
            lb(4 + (i-1)*n_coef_event) = pxSzT;
            ub(4 + (i-1)*n_coef_event) = ub_m1(i)-lb_m1(i);
            lb(5 + (i-1)*n_coef_event) = lb_m2(i);
            ub(5 + (i-1)*n_coef_event) = ub_m2(i);
            lb(6 + (i-1)*n_coef_event) = pxSzX;
            ub(6 + (i-1)*n_coef_event) = ub_m2(i)-lb_m2(i);
            lb(7 + (i-1)*n_coef_event) = pxSzX;
            ub(7 + (i-1)*n_coef_event) = ub_m2(i)-lb_m2(i);
            
             % Linear inequality constraints
            % sample mean of expModGaussian should be in image (m=mu+tau)
            A_ineq(2*i-1,6 + (i-1)*n_coef_event) = 1;
            A_ineq(2*i-1,8 + (i-1)*n_coef_event) = 1;
            A_ineq(2*i,6 + (i-1)*n_coef_event) = -1;
            A_ineq(2*i,8 + (i-1)*n_coef_event) = -1;
            b_ineq(2*i-1) = ub(6 + (i-1)*n_coef_event);
            b_ineq(2*i) = -lb(6 + (i-1)*n_coef_event);
            
        end
        lb(end) = min(imgPair(:));
        ub(end) = max(imgPair(:));
        
        % add baseline at the end
        p0 = [p0,bs];
        
        
    case 'tCaSpike_xEMG'
                       
        %coefficients: [t0, FM, tA, tI, FI, m_x, sd_x, expTauD_x]
        fun_event = @(p,x,y) ...
            (x < p(1)) .* 0 + ...
            (x >= p(1)) .* (0 + (p(2)*((-3*p(3)*(2 + (-1 + p(5))*p(3))*(4*p(3) - 7*p(4))*(p(3) - 3*p(4)))./(2*exp((2*(x-p(1)))/p(3)))+...
            (3*p(3)*(p(3) - 3*p(4))*(p(3) - 2*p(4)))./exp((x-p(1))/p(3)) +(3*p(3)*(5 + (-1 + p(5))*p(3))*(2*p(3) - 5*p(4))*(p(3) - p(4)))./(5*exp((5*(x-p(1)))/p(3)))-...
            (p(3)*(6 + (-1 + p(5))*p(3))*(p(3) - 2*p(4))*(p(3) - p(4)))./(6*exp((6.*(x-p(1)))/p(3))) -(6*p(4)^3*(1 + (-1 + p(5))*p(4)))./exp((x-p(1))/p(4))+...
            (18*p(4)^3*(p(3) + p(4) + (-1 + p(5))*p(3)*p(4)))./(exp((p(3)^(-1) + p(4)^(-1)).*(x-p(1)))*(p(3) + p(4))) -...
            (18*p(4)^3*(p(3) + (2 + (-1 + p(5))*p(3))*p(4)))./(exp((2/p(3) + p(4)^(-1)).*(x-p(1)))*(p(3) + 2*p(4)))+...
            (6*p(4)^3*(p(3) + (3 + (-1 + p(5))*p(3))*p(4)))./(exp((3/p(3) + p(4)^(-1)).*(x-p(1)))* (p(3) + 3*p(4))) -...
            (3*p(3)*(4 + (-1 + p(5))*p(3))*(5*p(3)^2 - 20*p(3)*p(4) +17*p(4)^2))./(4.*exp((4.*(x-p(1)))/p(3))) +...
            (p(3)*(p(3)^2*(57 + 10*(-1 + p(5))*p(3)) - 3*p(3)*(84 + 13*(-1 + p(5))*p(3))*p(4) + (249 + 29*(-1 + p(5))*p(3))*p(4)^2))./(3.*exp((3.*(x-p(1)))/p(3)))-...
            ((-1 + p(5))*(p(3) - 3*p(4))*(p(3) - 2*p(4))*(p(3) - p(4))*(37*p(3)^4 + 252*p(3)^3*p(4) + 605*p(3)^2*p(4)^2 + 660*p(3)*p(4)^3 + 360*p(4)^4))./...
            (60*(p(3) + p(4))*(p(3) + 2*p(4))*(p(3) + 3*p(4))) +  (6*(-1 + p(5))*p(3)^2*(p(3) - 3*p(4))*(p(3) - 2*p(4)).*...
            cosh((x-p(1))/p(3)))./exp((2.*(x-p(1)))/p(3))))/...
            ((p(3) - 3*p(4))*(p(3) - 2*p(4))*(p(3) - p(4)))).* ...
            exp(p(6)/p(8) + p(7)^2/(2*p(8)^2) - y./p(8)) .* ...
            cdf('Normal',y,p(6)+(p(7)^2/p(8)),p(7));
               
        n_coef_event = 8;
        coeff_n = {'t0','FM','tauR','tauD','FI',...
            'peakPos_x', 'sd_x', 'expTauD_x'};    
 
        p0 = zeros(1,N_events*n_coef_event);
        lb = zeros(1,N_events*n_coef_event+1);
        ub = inf(1,N_events*n_coef_event+1);
        A_ineq = zeros(2*N_events,N_events*n_coef_event+1);
        b_ineq = zeros(size(A_ineq,1),1);
    
        % parameters: [t0, FM, tA, tI, FI, m_x, sd_x, expTauD_x]
        for i=1:N_events
            
            p0(1 + (i-1)*n_coef_event) = prevFitParams.t0_est(i);
            p0(2 + (i-1)*n_coef_event) = prevFitParams.FM(i); % ...
                 % prevFitParams.maxOfEvents(i)*5;
            p0(3 + (i-1)*n_coef_event) = prevFitParams.tA(i);   %5;
            p0(4 + (i-1)*n_coef_event) = prevFitParams.tI(i);   %10;
            p0(5 + (i-1)*n_coef_event) = 1; %prevFitParams.FI(i);  
            p0(6 + (i-1)*n_coef_event) = prevFitParams.posOfEvents_x(i);
            p0(7 + (i-1)*n_coef_event) = prevFitParams.sd_x(i);
            p0(8 + (i-1)*n_coef_event) = prevFitParams.expTauD_x(i);
          
            % bounds for parameters
            lb(1 + (i-1)*n_coef_event) = min(sPP_pair(i),lb_m1(i) - pxSzT);
            try
                ub(1 + (i-1)*n_coef_event) = min(sPP_pair(i+1)-pxSzT,prevFitParams.posOfEvents_t(i));
            catch
                ub(1 + (i-1)*n_coef_event) = prevFitParams.posOfEvents_t(i);
            end
            
            lb(3 + (i-1)*n_coef_event) = pxSzT;
            ub(3 + (i-1)*n_coef_event) = ub_m1(i)-lb_m1(i);
            lb(4 + (i-1)*n_coef_event) = pxSzT;
            ub(4 + (i-1)*n_coef_event) = ub_m1(i)-lb_m1(i);
            lb(5 + (i-1)*n_coef_event) = 1;  
            ub(5 + (i-1)*n_coef_event) = 1;  
            lb(6 + (i-1)*n_coef_event) = lb_m2(i);
            ub(6 + (i-1)*n_coef_event) = ub_m2(i);
            lb(7 + (i-1)*n_coef_event) = pxSzX;
            ub(7 + (i-1)*n_coef_event) = ub_m2(i)-lb_m2(i);
            lb(8 + (i-1)*n_coef_event) = pxSzX;
            ub(8 + (i-1)*n_coef_event) = ub_m2(i)-lb_m2(i);
            
            % Linear inequality constraints
            % sample mean of expModGaussian should be in image (m=mu+tau)
            A_ineq(2*i-1,6 + (i-1)*n_coef_event) = 1;
            A_ineq(2*i-1,8 + (i-1)*n_coef_event) = 1;
            A_ineq(2*i,6 + (i-1)*n_coef_event) = -1;
            A_ineq(2*i,8 + (i-1)*n_coef_event) = -1;
            b_ineq(2*i-1) = ub(6 + (i-1)*n_coef_event);
            b_ineq(2*i) = -lb(6 + (i-1)*n_coef_event);
        end
        lb(end) = min(imgPair(:));
        ub(end) = max(imgPair(:));
        
        % add baseline at the end
        p0 = [p0,bs];
          
        % any(  (ub-lb)<0 )
        
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  keyboard
%         [F,F_events,F_individualEvents] = ...
%             sumEvents(p0,T,X,n_coef_event,fun_event,N_events);
%         
%         figure
%         mesh(T,X,imgPair,'LineStyle','-',...
%             'LineWidth',0.5,'FaceColor','none',...
%             'EdgeColor','k','EdgeAlpha',0.4)
%         hold on
%         mesh(T,X,F_individualEvents(:,:,1)+p0(end),'LineStyle','-',...
%             'LineWidth',0.5,'FaceColor','b',...
%             'EdgeColor','none','EdgeAlpha',0.1,'FaceAlpha',0.4)
%         
%         mesh(T,X,F_individualEvents(:,:,2)+p0(end),'LineStyle','-',...
%             'LineWidth',0.5,'FaceColor','g',...
%             'EdgeColor','none','EdgeAlpha',0.1,'FaceAlpha',0.4)
%         
%         mesh(T,X,F,'LineStyle','-',...
%             'LineWidth',0.5,'FaceColor','k',...
%             'EdgeColor','none','EdgeAlpha',0.1,'FaceAlpha',0.2)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
end

end


% sum N_events together using chosen function and proper initial
% parameters
function [F,F_events,F_individualEvents] = ...
    sumEvents(p,X,Y,n_coef_event,fun_event,N_events)

F_events = zeros(size(X));
F_individualEvents = zeros([size(X),N_events]);
for j=1:N_events
    
    F_event = fun_event(p(1+(j-1)*n_coef_event:n_coef_event+(j-1)*n_coef_event),X,Y);
    F_event(~isfinite(F_event))=0;
    
    F_events = F_events + F_event;
    F_individualEvents(:,:,j) = F_event;
    
end

F = F_events + p(end);

end


% fit function to use with fminsearch, fmincon; function returning
% scalar
function resSum = fitFunSumOfSquaers( ...
    p,X,Y,imgPair,n_coef_event,fun_event,N_events)
[F,~,~] = sumEvents(p,X,Y,n_coef_event,fun_event,N_events);
resSum = sum(sum((imgPair - F).^2));
end
% weighted
function resSum = fitFunSumOfSquaersW( ...
    p,X,Y,imgPair,n_coef_event,fun_event,N_events,W)

[F,~,~] = sumEvents(p,X,Y,n_coef_event,fun_event,N_events);
resSum = sum(sum(W.*(imgPair - F).^2));

end















