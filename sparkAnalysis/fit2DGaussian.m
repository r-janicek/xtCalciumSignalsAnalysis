function out = fit2DGaussian( ...
    X, Y, X_ups, Y_ups, imgE, imgE_m, peakPos, model)
keyboard
% X,Y,imgE = matrices
% imgE_m = mask of imgE, of detected event
% peakPos = pos of peak(row (y or spatial axis),column (x or time axis)) in D

% start and end of event in x direction (time)
sE = find(imgE_m(peakPos(1),:),1,'first');
eE = find(imgE_m(peakPos(1),:),1,'last');

pxSzT = mean(diff(X(1,:)));
pxSzX = mean(diff(Y(:,1)));

% guess parameters p0 for spark event
bs = prctile(imgE(~imgE_m),20); % baseline in D
A = max(imgE(imgE_m)); % amplitude

% axis for events, x and y
x_E = X(1,imgE_m(round(peakPos(1)),:));
y_E = Y(imgE_m(:,round(peakPos(2))),1);

% estimation of some p0
sd_x = std(x_E)/2; % std of gauss in x dimension
sd_y = std(y_E)/2; % std of gauss in y dimension

% profiles through max
x_prof_E = imgE(round(peakPos(1)),imgE_m(round(peakPos(1)),:)); 
x_prof_E = x_prof_E(:);
y_prof_E = imgE(imgE_m(:,round(peakPos(2))),round(peakPos(2))); 
y_prof_E = y_prof_E(:);

% function to get taus from profiles 
% rise
funExpRise = @(p,t,y) ((t>=p(1)).*((1-exp(-(t-p(1))./p(2))).*(p(3)) + p(4)) + ...
       (t<p(1)).*(p(4)))-y;
funExpRise_e = @(p,t) ((t>=p(1)).*((1-exp(-(t-p(1))./p(2))).*(p(3)) + p(4)) + ...
       (t<p(1)).*(p(4)));
% decay
funExpDecay = @(A,tD,B,x) A*exp(-x/tD) + B;

[~,mp_x_prof] = max(x_prof_E);
[~,mp_y_prof] = max(y_prof_E);

try
   % parts of profiles 
   % fit rise to get estimate of t0
   x_rise = x_E(1:mp_x_prof); x_rise = x_rise(:);
   x_prof_rise = x_prof_E(1:mp_x_prof); x_prof_rise = x_prof_rise(:);
   options = optimoptions('lsqnonlin','TolFun',10-9,'TolX',10-9,'MaxIter',1000,...
       'MaxFunEvals',3000,'Display','off');
   
   % t0,tR,A,y0
   x0 = [sE*pxSzT 3 A bs];
   pRise = lsqnonlin(@(p)funExpRise(p,x_rise,x_prof_rise),x0,...
       [zeros(1,3),min(x_prof_rise)],ones(1,4).*inf,options);
   t0_est = pRise(1);
   expTauR_x = pRise(2);
   
   if t0_est >= x_E(mp_x_prof-1)
       t0_est = sE*pxSzT;
   end
   if expTauR_x >= (x_E(mp_x_prof)-x_E(1))
       expTauR_x = 3;
   end
   
   % get estimate of decay taus
   x_decay = x_E(mp_x_prof:end); x_decay = x_decay(:);
   x_prof_decay = x_prof_E(mp_x_prof:end); x_prof_decay = x_prof_decay(:);
   y1_decay = y_E(mp_y_prof:end); y1_decay = y1_decay(:);
   y1_prof_decay = y_prof_E(mp_y_prof:end); y1_prof_decay = y1_prof_decay(:);
   y2_decay = y_E(1:mp_y_prof); y2_decay = y2_decay(:);
   y2_prof_decay = flipud( y_prof_E(1:mp_y_prof) ); y2_prof_decay = y2_prof_decay(:);
         
   coef_x = coeffvalues(fit(x_decay,x_prof_decay, funExpDecay, 'StartPoint',[A 10 bs],...
                'Lower',[0 0 0], ...
                'Upper',[inf (max(x_decay)-min(x_decay))*2 inf]));
            
   coef_y1 = coeffvalues(fit(y1_decay,y1_prof_decay, funExpDecay, 'StartPoint',[A 2 bs],...
                'Lower',[0 0 0], 'Upper',[inf inf inf]));  
            
   coef_y2 = coeffvalues(fit(y2_decay,y2_prof_decay, funExpDecay, 'StartPoint',[A 2 bs],...
                'Lower',[0 0 0], 'Upper',[inf inf inf]));
                 
   expTauD_x = coef_x(2);  % tau decay in x
   expTauD_y = mean( [coef_y1(2), coef_y2(2)] );  % tau decay in y
   
   if expTauD_x >= (x_E(end)-x_E(mp_x_prof))
       expTauD_x = (x_E(end)-x_E(mp_x_prof));
   end
   if expTauD_y >= (y_E(end)-y_E(1))*0.68
       expTauD_y = (y_E(end)-y_E(1))*0.68;
   end
   
   
catch
    t0_est = sE*pxSzT;
    expTauR_x = 3;
    expTauD_x = 10;
    expTauD_y = 1.5;
end

% position of gauss peak in EMG in x direction
m1 = peakPos(2)*pxSzT-expTauD_x/4; if m1<pxSzT, m1=pxSzT; end
% position of gauss peak in EMG in y direction
m2 = peakPos(1)*pxSzX-expTauD_y/2; if m2<pxSzX, m2=pxSzX; end 

% bounds for positions parameters
lb_m1 = sE*pxSzT-pxSzT; if lb_m1<0, lb_m1=0; end
ub_m1 = eE*pxSzT;
lb_m2 = min(min(Y(imgE_m)));
ub_m2 = max(max(Y(imgE_m)));

switch model
    
%     case '2D_G'
%         % 2D gaussian
%         % coefficients: [A, mu_x, mu_y, sd_x, sd_y F0]
%         p0 = [A peakPos(2)*pxSzT peakPos(1)*pxSzX sd_x sd_y bs];
%         
%         fitFun = @(p,x,y) p(1).*exp(-((x-p(2)).^2./(2*p(4).^2))) .* exp(-((y-p(3)).^2./(2*p(5).^2))) + p(6);
%         
% add lower and upper bonds for positions in x and y coordinates
%         lb = zeros(size(p0));
%         ub = inf(size(p0));
%         lb(2) = lb_m1;
%         ub(2) = ub_m1;
%         lb(3) = lb_m2;
%         ub(3) = ub_m2;
%         lb(4) = pxSzT;
%         ub(4) = ub_m1-lb_m1;
%         lb(5) = pxSzX;
%         ub(5) = ub_m2-lb_m2;
%         lb(end) = min(imgE(:));
         
              
    case '2D_EMGt'
        % 2D exponentially modified gaussian in time direction only
        % coefficients: [A, mu_x, sd_x, expTauD, mu_y, sd_y, F0]       
        p0 = [A m1 sd_x expTauD_x peakPos(1)*pxSzX sd_y bs];  
        
        fitFun = @(p,x,y) p(1) .* exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - x./p(4)) .* ...
            cdf('Normal',x,p(2)+(p(3)^2/p(4)),p(3)) .* ...
            exp(-((y-p(5)).^2./(2*p(6).^2))) + p(7);
        
        % add lower and upper bonds for positions in x and y coordinates
        lb = zeros(size(p0));
        ub = inf(size(p0));
        lb(2) = lb_m1;
        ub(2) = ub_m1;
        lb(3) = pxSzT;
        ub(3) = ub_m1-lb_m1;
        lb(4) = pxSzT;
        ub(4) = ub_m1-lb_m1;
        lb(5) = lb_m2;
        ub(5) = ub_m2;
        lb(6) = pxSzX;
        ub(6) = ub_m2-lb_m2;
        lb(end) = min(imgE(:));
        
        % Linear inequality constraints
        % sample mean of expModGaussian should be in image (m=mu+tau)    
        A_ineq = zeros(2,numel(p0));
        b_ineq = zeros(size(A_ineq,1),1);
        A_ineq(1,[2,4]) = 1;
        A_ineq(2,[2,4]) = -1;

        b_ineq(1) = ub(2);
        b_ineq(2) = -lb(2);
       
        
                
    case '2D_EMGxt'
        % 2D exponentially modified gaussian in each direction
        % coefficients: [A, mu_x, sd_x, expTauD_x, mu_y, sd_y, expTauD_y, F0]
        p0 = [sd_x*sd_y*A m1 sd_x expTauD_x m2 sd_y expTauD_y bs];
        
        fitFun = @(p,x,y) p(1) .* exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - x./p(4)) .* ...
            cdf('Normal',x,p(2)+(p(3)^2/p(4)),p(3)) .* ...
            exp(p(5)/p(6) + p(6)^2/(2*p(7)^2) - y./p(7)) .* ...
            cdf('Normal',y,p(5)+(p(6)^2/p(7)),p(6)) + p(8);
        
        % add lower and upper bonds for positions in x and y coordinates
        lb = zeros(size(p0));
        ub = inf(size(p0));
        lb(2) = lb_m1;
        ub(2) = ub_m1;
        lb(3) = pxSzT;
        ub(3) = ub_m1-lb_m1;
        lb(4) = pxSzT;
        ub(4) = ub_m1-lb_m1;
        lb(5) = lb_m2;
        ub(5) = ub_m2;
        lb(6) = pxSzX;
        ub(6) = ub_m2-lb_m2;
        lb(7) = pxSzX;
        ub(7) = ub_m2-lb_m2;
        lb(end) = min(imgE(:));
        
        % [A, mu_x, sd_x, expTauD_x, mu_y, sd_y, expTauD_y, F0]
        % Linear inequality constraints
        % sample mean of expModGaussian should be in image (m=mu+tau)    
        A_ineq = zeros(4,numel(p0));
        b_ineq = zeros(size(A_ineq,1),1);
        A_ineq(1,[2,4]) = 1;
        A_ineq(2,[2,4]) = -1;
        A_ineq(3,[5,7]) = 1;
        A_ineq(4,[5,7]) = -1;
       
        b_ineq(1) = ub(2);
        b_ineq(2) = -lb(2);
        b_ineq(3) = ub(5);
        b_ineq(4) = -lb(5);
        
        
    case 'tCaSpike_xEMG'

        %coefficients: [t0, FM, tA, tI, FI, mu_x, sd_x, expTauD_x, F0]
        p0 = [t0_est A expTauR_x expTauD_x 1 m2 sd_y expTauD_y bs];
        fitFun = @(p,x,y) ...
            (x < p(1)) .* p(9) + ...
            (x >= p(1)) .* (p(9) + ((p(2)*((-3*p(3)*(2 + (-1 + p(5))*p(3))*(4*p(3) - 7*p(4))*(p(3) - 3*p(4)))./(2*exp((2*(x-p(1)))/p(3)))+...
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
            cdf('Normal',y,p(6)+(p(7)^2/p(8)),p(7)));
 
        % add lower and upper bonds for positions in x and y coordinates
        % coefficients: [t0, FM, tA, tI, FI, m_x, sd_x, expTauD_x, F0]
        lb = zeros(size(p0));
        ub = inf(size(p0));
        lb(1) = lb_m1 - 2*pxSzT;
        ub(1) = peakPos(2)*pxSzT;
        lb(3) = pxSzT;
        ub(3) = ub_m1-lb_m1;
        lb(4) = pxSzT;
        ub(4) = ub_m1-lb_m1;
        lb(5) = 1;
        ub(5) = 1;
        lb(6) = lb_m2;
        ub(6) = ub_m2;
        lb(7) = pxSzX;
        ub(7) = ub_m2-lb_m2;
        lb(8) = pxSzX;
        ub(8) = ub_m2-lb_m2;
        lb(end) = min(imgE(:));
        
        % Linear inequality constraints
        % sample mean of expModGaussian should be in image (m=mu+tau)    
        A_ineq = zeros(2,numel(p0));
        b_ineq = zeros(size(A_ineq,1),1);
        A_ineq(1,[6,8]) = 1;
        A_ineq(2,[6,8]) = -1;
       
        b_ineq(1) = ub(6);
        b_ineq(2) = -lb(6);
               
%         figure
%         plot3(X,Y,imgE,'LineStyle','none','Marker','.','Color',[0.1 0.1 0.1])
%         hold on
%         surf(X,Y,fitFun(p0,X,Y),'FaceAlpha',0.6,'EdgeColor','none','FaceColor','interp')
      
end
 
% first fit with only event defined by mask, zero weights for not event
% points
wE = zeros(size(imgE));
wE(imgE_m) = ( 1 + imgE(imgE_m)./A ).^4;
%wE = ones(size(imgE));
% sum of squares, scalar
fitFunSum = @(p,X,Y,D,wE) sum( sum( wE.*(D - fitFun(p,X,Y)).^4 ) );

% options for fmincon
opt = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
    'MaxIter',1000,'MaxFunEvals',3000,'UseParallel',false);

try
    coef2DGauss = fmincon(@(p)fitFunSum(p,X,Y,imgE,wE), ... 
        p0,A_ineq,b_ineq,[],[],lb,ub,[],opt);
catch
    coef2DGauss = p0;
end

% fit baseline before the start of first event
try
    wE = ones(size(imgE));
    % start of pair of events
    ind_E_s = find(sum(imgE_m,1)>0,1,'first');
    wE(:,ind_E_s:end) = 0;

    lb_bs = [coef2DGauss(1:end-1),min(imgE(:))];
    ub_bs = [coef2DGauss(1:end-1),max(imgE(:))];
    coef2DGauss = fmincon(@(p)fitFunSum(p,X,Y,imgE,wE),...
        coef2DGauss,A_ineq,b_ineq,[],[],lb_bs,ub_bs,[],opt);
catch
 
end

% refit with higher weights set to event area and zero otherwise, also with
% locked baseline
wE = zeros(size(imgE));
wE(imgE_m) = ( 1 + imgE(imgE_m)./A ).^2;
lb(end) = coef2DGauss(end);
ub(end) = coef2DGauss(end);

try
    coef2DGauss = fmincon(@(p)fitFunSum(p,X,Y,imgE,wE), ...
        coef2DGauss,A_ineq,b_ineq,[],[],lb,ub,[],opt);
catch 
    try
        coef2DGauss = fmincon(@(p)fitFunSum(p,X,Y,imgE,wE), ...
            zeros(size(p0)),[],[],[],[],-inf(size(p0)),inf(size(p0)),[],opt);
    catch 
        coef2DGauss = p0;
    end
end

D_fit = fitFun(coef2DGauss,X,Y);
D_fit_ups = fitFun(coef2DGauss,X_ups,Y_ups);

% create output structure
out = struct('coef2DGauss',coef2DGauss,...
             'X',X,'Y',Y,'dataEvent',imgE,...
             'dataEventFit',D_fit,...
             'X_ups',X_ups,'Y_ups',Y_ups,...
             'dataEventFit_ups',D_fit_ups);

         
%       keyboard    
% figure
% mesh(X,Y,imgE,'LineStyle','-',...
%     'LineWidth',0.5,'FaceColor','none',...
%     'EdgeColor','k','EdgeAlpha',0.4)
% hold on
% surf(X,Y,fitFun(coef2DGauss,X,Y),'FaceAlpha',0.6,'EdgeColor','none','FaceColor','interp')
% %                   
end



