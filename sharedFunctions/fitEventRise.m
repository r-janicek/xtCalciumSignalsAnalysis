function [h_line,coef,sp_fit,startOfSpark] = fitEventRise(...
    pxSzT,x_t,prof_t,prof_s,pks,locs,ax_prof,coefPrevFit,tol,iter,...
    smooth_span,bs_crit,posOfSelPoints)

options = optimoptions('lsqnonlin','TolFun',tol,'TolX',tol,'MaxIter',iter,...
    'MaxFunEvals',3*iter,'Display','off');

fun_e = @(x,t) ((t>=x(1)).*((1-exp(-(t-x(1))./x(2))).*(x(3)) + x(4)) + ...
        (t<x(1)).*(x(4)));
    
fun = @(x,t,ys) ((t>=x(1)).*((1-exp(-(t-x(1))./x(2))).*(x(3)) + x(4)) + ...
         (t<x(1)).*(x(4)))-ys;

if ~isempty(posOfSelPoints)
    pks = posOfSelPoints.peakFitVal;
    locs = x_t(posOfSelPoints.peakFitPos);
end

h_line = zeros(length(pks),1);
coef = zeros(length(pks),4); 
sp_fit = {zeros(length(pks),1)};
startOfSpark = zeros(length(pks),1); 
prof_s = prof_s(:);

maxDurOfBaseline = round(500/pxSzT); % maximum duration of baseline in points

% smooth profile for difference analysis, loess with defined duration in ms  
%prof_s = smooth(prof_t,3);
% n_pts = round(smooth_span/pxSzT);
% prof_s = smooth(prof_t,n_pts/length(prof_t),'loess');

bs_crit = round(bs_crit);

%remove baseline, only events bigger than specified percentile stay
percentl = prctile(prof_s,[25 50 100-bs_crit]);
%iqr = percentl(3)- percentl(1);
trsh = percentl(3); % threshold

prof_s(prof_s < trsh) = trsh;

% location of peaks in points
[~, locsPx] = min(abs(x_t'-locs)); % using expansion of minus to matrix, not working with old versions of matlab

prof_s_d = diff(prof_s);
prof_s_d = [0;prof_s_d];

% %find peaks of derivation of profile
% [vals_d_pos,locs_d_pos] = findpeaks(prof_s_d,x_t);

% look for positive peaks of derivation 20 ms around peak of spark
% d = 100;

if ~isempty(pks)
    
    for i=1:length(pks)
        
        %first try to find start position of baseline of spark with derivation
        try
            % take part of profile derivation
            testProfDer = zeros(size(prof_s_d));
            if i==1
                testProfDer(1:locsPx(i)) = prof_s_d(1:locsPx(i));
            else
                testProfDer(locsPx(i-1):locsPx(i)) = prof_s_d(locsPx(i-1):locsPx(i));
            end
            
            [~,pm] = max(testProfDer);
            
            % flip it
            testProfDer = flipud(testProfDer);
            [~,pmFl] = max(testProfDer);
            
            % find first negative derivation
            p_d_neg = find(testProfDer(pmFl:end)<0,1,'first');
            
            % corrected position for not flipped max
            p_d_neg = pm-p_d_neg;
            
            if isempty(p_d_neg)
                p_d_neg = 1;
            end
            
            pos_s = locsPx(i) - maxDurOfBaseline;
            
            if pos_s <= 0
                pos_s = 1;
            end
            
            if ~isempty(p_d_neg) && (p_d_neg > pos_s) && (p_d_neg<locsPx(i))
                pos_s = p_d_neg+1;
            end
            
           
            % if there is problem, use old approach
        catch
            
            if i==1
                
                [~,pos_m1] = min(abs(x_t-locs(i)));
                pos_s = pos_m1 - maxDurOfBaseline;
                
                if pos_s <= 0
                    pos_s = 1;
                end
                
            else
                
                [~,pos_m1] = min(abs(x_t-locs(i-1)));
                [~,pos_m2] = min(abs(x_t-locs(i)));
                % baseline of spark start in the middle between two  sparks
                dist12 = round((pos_m2 - pos_m1)*0.5);
                
                if dist12 > maxDurOfBaseline
                    dist12 = maxDurOfBaseline;
                    
                elseif dist12 < 5
                    dist12=5;
                end
                
                pos_s = pos_m2 - dist12;
                
            end
            
        end
        
        [~,pos_e] = min(abs(x_t-locs(i)));
        
        if ~isempty(posOfSelPoints)
            pos_s = posOfSelPoints.bsFitS;
            pos_e = posOfSelPoints.peakFitPos;
        end
        
        ys = prof_t(pos_s:pos_e);
        t = x_t(pos_s:pos_e);
         
        if length(ys)<size(coef,2)
            ys = prof_t(pos_s:pos_e + (size(coef,2)-length(ys)));
            t = x_t(pos_s:pos_e + (size(coef,2)-length(t)));
        end
        
        t = t(:);
        ys = ys(:);
        
        % estimate initial fit parameters
        % t0,tR,A,y0
        if ~isempty(coefPrevFit)
            x0 = coefPrevFit(i,:);
        else
            A = pks(i);
            y0 = prctile(ys,10);
            
            % estimate t0
            lastNegBs = find(ys-y0<0,1,'last');
            if lastNegBs >= numel(ys)
                lastNegBs = lastNegBs-1;
            end
            if lastNegBs < 1
                lastNegBs = 1;
            end
            if isempty(lastNegBs)
                lastNegBs = 1;
            end
            [~,maxInYs] = max(ys(lastNegBs:end));
           
            % find 25 and 75% from max in interval (lastNegBs:end)
            ysRisePart = ys(lastNegBs:lastNegBs+maxInYs-1);
            [~,p_75] = min( abs(ysRisePart-max(ysRisePart).*0.75) );
            [~,p_25] = min( abs(ysRisePart-max(ysRisePart).*0.25) ); 
            if p_75 == p_25
                p_25 = p_75-1;
            end
            p_25 = p_25+lastNegBs;
            p_75 = p_75+lastNegBs;
       
            v_75 = ys(p_75);
            v_25 = ys(p_25);
      
            % construct line two-point form and get x at y=0
            t0_ind_est = round( p_25 + ( (0-v_25)*(p_75-p_25) ) / (v_75-v_25) );
            tR_est = t(p_75)-t(p_25);
            if t0_ind_est<1, t0_ind_est=1; end
            if t0_ind_est>numel(t), t0_ind_est = numel(t); end
            if tR_est<1, tR_est=5; end
            
            % t0,tR,A,y0
            x0 = [t(t0_ind_est) tR_est A y0];
        end

        if ~isempty(posOfSelPoints)
            % t0,tR,A,y0
            x0(4) = mean(ys(1:posOfSelPoints.bsFitE-pos_s+1));
            % lower and upper bounds
            lb_fit = [zeros(1,3),x0(4)];
            ub_fit = [ones(1,3).*inf,x0(4)];
        else
            % lower and upper bounds
            lb_fit = [zeros(1,3),min(ys)];
            ub_fit = ones(1,4).*inf;
        end
        
        % fit
        try
            x = lsqnonlin(@(x)fun(x,t,ys),x0,lb_fit,ub_fit,options);
        catch
            x = x0;
        end
        
        coef(i,:) = x;
        sp_fit(i,1) = {[t,fun_e(x,t)]};
        startOfSpark(i,1) = pos_s; % in pixels
        
        if ~isempty(ax_prof)
            h_line(i,1) = line(t,fun_e(x,t),'Parent',ax_prof,'Color',[0.2 1 0.2],'LineStyle','-','LineWidth',2,'Tag','fitLineRise');
        end
        
    end
end

end


