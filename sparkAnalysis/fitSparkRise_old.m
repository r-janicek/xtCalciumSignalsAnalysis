function [h_line,detectedEventsMask,coef,sp_fit,startOfSpark,endOfSpark] = ...
    fitSparkRise(pxSz_t,x_t,prof_t,pks,locs,ax_prof,coefPrevFit,tol,iter,...
    smooth_span,bs_crit,sSpPrev,eSpPrev)

% locs in time units
options = optimoptions('lsqnonlin','TolFun',tol,'TolX',tol,'MaxIter',iter,...
    'MaxFunEvals',3*iter,'Display','off');

fun_e = @(x,t) ((t>=x(1)).*((1-exp(-(t-x(1))./x(2))).*(x(3)) + x(4)) + ...
        (t<x(1)).*(x(4)));
    
fun = @(x,t,ys) ((t>=x(1)).*((1-exp(-(t-x(1))./x(2))).*(x(3)) + x(4)) + ...
         (t<x(1)).*(x(4)))-ys;

h_line = zeros(length(pks),1);
coef = zeros(length(pks),4); 
sp_fit = {zeros(length(pks),1)};
startOfSpark = zeros(length(pks),1); 
endOfSpark = zeros(length(pks),1);
detectedEventsMask = false(numel(prof_t),1);

maxDurOfBaseline = round(50/pxSz_t); % maximum duration of baseline in points

% use previous starts and ends of sparks
if isempty(sSpPrev) && isempty(eSpPrev)
    
    % smooth profile for difference analysis, loess with defined duration in ms
    %prof_s = smooth(prof_t,3);
    n_pts = round(smooth_span/pxSz_t);
    prof_s = smooth(prof_t,n_pts/length(prof_t),'loess');
    
    bs_crit = round(bs_crit);
    
    %remove baseline, only events bigger than specified percentile stay
    percentl = prctile(prof_s,[25 50 bs_crit]);
    %iqr = percentl(3)- percentl(1);
    
    prof_s(prof_s < percentl(3))=percentl(3);
    
%       figure
%       plot(prof_s)
    
    [v_s,l_s] = findpeaks(prof_s);
    % indx = find(v_s<std(prof_t)+percentl(1));
    %
    % v_s(indx)=[];
    % l_s(indx)=[];
    
    bb = percentl(3);
    
    if ~isempty(v_s)
        for i=1:length(v_s)
            
            if i==1
                
                if i==length(v_s)
                    
                    pp1 = prof_s(1:l_s(i));
                    %calculate gradient from right to left
                    f_pp = flipud(pp1);
                    g_pp_r = gradient(f_pp);
                    [st_e,~] =find((g_pp_r>0)|(g_pp_r==0),1,'first');
                    st_e = length(pp1)-(st_e-2);
                    
                    pp1(1:st_e)=bb;
                    
                    %calculate gradient from left to right
                    pp2 = prof_s(l_s(i):end);
                    g_pp_l2 = gradient(pp2);
                    [end_e2,~] =find((g_pp_l2>0)|(g_pp_l2==0),1,'first');
                    end_e2 = end_e2-1;
                    
                    pp2(end_e2:end)=bb;
                    
                    pp = [pp1;pp2(2:end)];
                    
                else
                    
                    pp = prof_s(1:l_s(1));
                    f_pp = flipud(pp);
                    g_pp = gradient(f_pp);
                    [st_e,~] = find((g_pp>0)|(g_pp==0),1,'first');
                    st_e = length(pp)-(st_e-2);
                    pp(1:st_e)=bb;
                    
                end
                
            elseif i>1 && i==length(v_s)
                
                pp1 = prof_s(l_s(i-1):l_s(i));
                %calculate gradient from left to right
                g_pp_l = gradient(pp1);
                [end_e,~] =find((g_pp_l>0)|(g_pp_l==0),1,'first');
                end_e = end_e-1;
                
                %calculate gradient from right to left
                f_pp = flipud(pp1);
                g_pp_r = gradient(f_pp);
                [st_e,~] =find((g_pp_r>0)|(g_pp_r==0),1,'first');
                st_e = length(pp1)-(st_e-2);
                
                if  st_e > end_e
                    pp1(end_e:st_e)=bb;
                end
                
                pp2 = prof_s(l_s(i):end);
                g_pp_l2 = gradient(pp2);
                [end_e2,~] =find((g_pp_l2>0)|(g_pp_l2==0),1,'first');
                end_e2 = end_e2-1;
                
                pp2(end_e2:end)=bb;
                
                pp = [pp1;pp2(2:end)];
                
            else
                
                pp = prof_s(l_s(i-1):l_s(i));
                %calculate gradient from left to right
                g_pp_l = gradient(pp);
                [end_e,~] =find((g_pp_l>0)|(g_pp_l==0),1,'first');
                end_e = end_e-1;
                
                %calculate gradient from right to left
                f_pp = flipud(pp);
                g_pp_r = gradient(f_pp);
                [st_e,~] =find((g_pp_r>0)|(g_pp_r==0),1,'first');
                st_e = length(pp)-(st_e-2);
                
                if  st_e > end_e
                    pp(end_e:st_e)=bb;
                end
                
            end
            
            if exist('prof_s_n','var')
                prof_s_n = [prof_s_n;pp(2:end)];
            else
                prof_s_n = pp;
            end
            
        end
    end
  
    % calculate derivative
    prof_s_d = diff(prof_s_n);
    
    % find peaks of derivation of profile
    [vals_d_pos,locs_d_pos] = findpeaks(prof_s_d,x_t(2:end));
    
%     figure
%     plot(prof_s_d)
    
    % look for positive peaks of derivation 20 ms around peak of spark
    d = 20;
    
end


if ~isempty(pks)
   
    for i=1:length(pks)
            
        if isempty(sSpPrev) && isempty(eSpPrev)
            
            % first try to find start position of baseline of spark with derivation
            try
                
                if i==1
                    
                    % find position of max of 1.derivative which correspond to chosen peak
                    p_d_pos1 = find( (locs_d_pos < locs(i)) & (locs_d_pos > locs(i)-d) );
                    
                    p_d_v = vals_d_pos(p_d_pos1);
                    p_d_pos1 = locs_d_pos(p_d_pos1);
                    
                    [~,indx_pos] = sort(p_d_v,'descend');
                    p_d_pos1 = p_d_pos1(indx_pos);
                    p_d_pos = p_d_pos1(1);
                    
                    [~,p_d_pos] = min(abs(x_t-p_d_pos));
                    
                    % find last negative derivative/end of previous spark
                    p_d_neg = find(prof_s_d(1:p_d_pos) < 0,1,'last');
                    
                    % find position of point 50 ms before peak 
                    % and take that if there is no negative derivative before
                    [~,pos_peak] = min(abs(x_t-locs(i))); % in points
                    
                    pos_s = pos_peak - maxDurOfBaseline;
                    
                    if pos_s <= 0
                        pos_s = 1;
                    end
                    
                    if ~isempty(p_d_neg) && (p_d_neg > pos_s) && (p_d_neg<pos_peak)
                        pos_s = p_d_neg+1;
                    end
                    
                    % find end of event, look in derivation of profile
                    % find first negative derivative and look for first
                    % positive derivative
                    posNeg1 = find(prof_s_d(p_d_pos:end)<0,1,'first')+p_d_pos-1;
                    posPos1 = find(prof_s_d(posNeg1:end)>=0,1,'first')+posNeg1-1;
                    
                    pos_e = pos_peak + 2*maxDurOfBaseline;
                    
                    if pos_e >= length(prof_t)
                        pos_e = length(prof_t);
                    end
                    
                    if ~isempty(posPos1) && (posPos1 < pos_e) && (posPos1>pos_peak)
                        pos_e = posPos1;
                    end
                    
                    
                else
                    
                    % start of event
                    if (locs(i)-d) > locs(i-1)
                        p_d_pos1 = find( (locs_d_pos < locs(i)) & (locs_d_pos > locs(i)-d) );
                    else
                        p_d_pos1 = find( (locs_d_pos > locs(i-1)) & (locs_d_pos < locs(i)) );
                    end
                    
                    p_d_v = vals_d_pos(p_d_pos1);
                    p_d_pos1 = locs_d_pos(p_d_pos1);
                    
                    [~,indx_pos] = sort(p_d_v,'descend');
                    p_d_pos1 = p_d_pos1(indx_pos);
                    p_d_pos = p_d_pos1(1);
                    
                    [~,p_d_pos] = min(abs(x_t-p_d_pos));
                    
                    p_d_neg = find(prof_s_d(1:p_d_pos) < 0,1,'last');
                    
                    [~,pos_peak] = min(abs(x_t-locs(i)));  % in points
                    
                    pos_s = pos_peak - maxDurOfBaseline;
                    
                    if ~isempty(p_d_neg) && (p_d_neg > pos_s) && (p_d_neg<pos_peak)
                        pos_s = p_d_neg+1;
                    end
                    
                    % end of event
                    posNeg1 = find(prof_s_d(p_d_pos:end)<0,1,'first')+p_d_pos-1;
                    posPos1 = find(prof_s_d(posNeg1:end)>=0,1,'first')+posNeg1-1;
                                        
                    pos_e = pos_peak + 2*maxDurOfBaseline;
                    
                    if pos_e >= length(prof_t)
                        pos_e = length(prof_t);
                    end
                    
                    if ~isempty(posPos1) && (posPos1 < pos_e) && (posPos1>pos_peak)
                        pos_e = posPos1;
                    end
                    
                end
                
                % if there is problem, use old approach
            catch
                
                if i==1
                    
                    [~,pos_m1] = min(abs(x_t-locs(i)));
                    pos_s = pos_m1 - maxDurOfBaseline;
                    
                    if pos_s <= 0
                        pos_s = 1;
                    end
                    
                    pos_e = pos_m1 + 2*maxDurOfBaseline;
                    
                    if pos_e >= length(prof_t)
                        pos_e = length(prof_t);
                    end
                    
                else
                    
                    [~,pos_m1] = min(abs(x_t-locs(i-1)));
                    [~,pos_m2] = min(abs(x_t-locs(i)));
                    
                    % baseline of event start in the middle between two events
                    dist12 = round((pos_m2 - pos_m1)*0.5);
                    
                    if dist12 > maxDurOfBaseline
                        dist12 = maxDurOfBaseline;
                        
                    elseif dist12 < 5
                        dist12=5;
                    end
                    
                    pos_s = pos_m2 - dist12;
                    
                    % end of event
                    % baseline of event ends in the middle between two events
                    try
                        [~,pos_m3] = min(abs(x_t-locs(i+1)));
                        dist23 = round((pos_m3 - pos_m2)*0.5);
                        
                        if dist23 > 2*maxDurOfBaseline
                            dist23 = 2*maxDurOfBaseline;
                        end
                        
                        pos_e = pos_m3 - dist23;
                        
                    catch
                        pos_e = pos_m2 + 2*maxDurOfBaseline;
                        
                        if pos_e >= length(prof_t)
                            pos_e = length(prof_t);
                        end
                        
                    end
                    
                end
                
            end
        
            
        else
            % take old
            pos_s = sSpPrev(i);
            pos_e = eSpPrev(i);
            
        end
        
        
        % get position of peak 
        [~,pos_p] = min(abs(x_t-locs(i)));
        
        ys = prof_t(pos_s:pos_p);
        t = x_t(pos_s:pos_p);

        if length(ys)<size(coef,2)
            
            ys = prof_t(pos_s:pos_p + (size(coef,2)-length(ys)));
            t = x_t(pos_s:pos_p + (size(coef,2)-length(t)));
            
        end
        
        t = t(:);
        ys = ys(:);
        
                       
        % t0,tR,A,y0
        if ~isempty(coefPrevFit)
            x0 = coefPrevFit(i,:);
        else
            % estimation of initial fit values
            A = pks(i);
            
            t0FirstEst = locs(i)-10;
            [~,t0FirstEstPx] = min( abs( t-t0FirstEst ) );
            if t0FirstEstPx <= 0
                t0FirstEstPx = 3; % first five points
            end
            y0 = mean(ys(1:t0FirstEstPx));
            try
                % estimate t0
                % find 25 and 75% from max of rise part
                [~,p_85] = min( abs( ys-(y0 + max(ys).*0.75) ) );
                [~,p_35] = min( abs( ys-(y0 + max(ys).*0.25) ) );
                if p_85 == p_35
                    p_35 = p_85-1;
                end
                
                v_85 = ys(p_85);
                v_35 = ys(p_35);
                
                % construct line two-point form and get x at y=y0
                t0_ind_est = round( p_35 + ( (y0-v_35)*(p_85-p_35) ) / (v_85-v_35) );
                tR_est = t(p_85)-t(p_35);
                if t0_ind_est<1, t0_ind_est = 1; end
                % t0,tR,A,y0
                x0 = [t(t0_ind_est) tR_est A-y0 y0];
            catch
                % t0,tR,A,y0
                x0 = [locs(i)-10 5 A y0];
            end
        end

        try
            x = lsqnonlin(@(x)fun(x,t,ys),x0,[zeros(1,3),min(ys)],ones(1,4).*inf,options);       
        catch
            try
                sumSqrs = @(x,t,ys) sum(fun(x,t,ys).^2);
                x = fmincon(@(x)sumSqrs(x,t,ys),x0,[],[],[],[],...
                    [zeros(1,3),min(ys)],ones(1,4).*inf); 
            catch 
               keyboard 
            end
        end
         
%         figure
%          plot(t,ys,'ok')
%          hold on 
%          plot(t,fun_e(x0,t),'ob')
%           
%          plot(t,fun_e(x,t),'or')
          
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












