function eventParams = calcEventParams(Img_data,S_event,sparkStart,prevFitCoef,posOfROI,main_fig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% dorobit este pre sparky iba
% aby sa dala rozsirit oblast pre spark, tak aby sa nebral do uvahy
% susedny spark ak je tam nejaky; mozno skusit fciu ismember alebo nieco podobne, mozno rekurzia
% alebo ak je to moc komplikovane len skusit bez rozsirovania
%

if isempty(sparkStart)
    sparkStart = zeros(1,numel(S_event));
end

%get type of analysis
analysis_type = getappdata(main_fig,'analysisType');
pxSize_t = getappdata(main_fig,'pxSize_t');
pxSz = getappdata(main_fig,'pxSize_x');
%smooth_span = str2double(get(getappdata(main_fig,'h_smooth_edit'),'String'));
%n_pts = round(smooth_span/pxSize_t);
scrRes = get(0,'ScreenSize');
% to calculate parameters take x and t average profiles from event
% diameter around event centre - from the GUI for time profile
% 5 ms around max position (after smoothing with csaps) for x profile
h_d = str2double(get(getappdata(main_fig,'h_edit_averageWidth'),'String'));

n_px_t = ceil(h_d/pxSz);
if mod(n_px_t,2)==0
    n_px_t = n_px_t - 1;
end

n_px_x = round(5/pxSize_t);
if mod(n_px_x,2)==0
    n_px_x = n_px_x - 1;
end

% pozret to este, aby som vzdy dostal parametre aj ked fit zlyha, plus
% pouzit parametre z predchadzajuceho fitovania

for i = 1:length(S_event)
    %keyboard
    % try to expand area around spark, in time use start of spark from profile
    % in spatial axis 2 um; to get baseline
    rows_e = S_event(i).SubarrayIdx{1};
    cols_e = S_event(i).SubarrayIdx{2};
    
    
    %     keyboard
    %    %%%%%%%%%%%%%%%%%
    %    % try to expand spark area based on ismember and pixels position of
    %    % others events in image data
    %     testEventPixId = S_event(i).SubarrayIdx;
    %     othersE = S_event;
    %     othersE(i)=[];
    %
    %     othersEventsPixId = vertcat(othersE(:).PixelIdxList);
    %
    %
    %     ismember(testEventPixId,othersEventsPixId)
    %
    %    %%%%%%%%%%%%%%%%%
    
    % do not do extension for sparks detection
    if ~strcmp(analysis_type,'spark detection')
        
        ext_x = round(2/pxSz);   % extension in x, in pixels
        
        ext_t_s = cols_e(1) - sparkStart(i); % extension in t, in pixels
        
        % set the extension of end of spark in pixels, fixed to 100 ms if more
        % or last one
        if i==length(S_event)
            ext_t_e = round(100/pxSize_t);
        else
            ext_t_e = sparkStart(i+1)-cols_e(end);
            
            if ext_t_e*pxSize_t > 100
                
                ext_t_e = round(100/pxSize_t);
                
            end
        end
        
        %check if it is OK in terms of size of Img_data
        if (size(Img_data,1)<rows_e(end)+ext_x) || (rows_e(1)-ext_x <= 0)
            
            ext_x = min(size(Img_data,1)-rows_e(end)-1, rows_e(1)-1);
            
        end
        
        if (cols_e(1)-ext_t_s) <= 0
            
            ext_t_s = cols_e(1)-1;
            
        end
        
        if size(Img_data,2) < (cols_e(end)+ext_t_e)
            
            ext_t_e = size(Img_data,2)-cols_e(end)-1;
            
        end
        
        %expand spark area
        rows_e = [(rows_e(1)-ext_x:rows_e(1)-1),rows_e,(rows_e(end)+1:rows_e(end)+ext_x)];
        cols_e = [(cols_e(1)-ext_t_s:cols_e(1)-1),cols_e,(cols_e(end)+1:cols_e(end)+ext_t_e)];
    end
    
    D = Img_data(rows_e,cols_e);
    
    %     [X,Y] = meshgrid((1:1:size(D,2)), (1:1:size(D,1)));
    %     [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(X,Y,D)
    %     figure
    %     scatter3(X(:),Y(:),D(:))
    %     hold on
    %     surf(X,Y,zfit)
    
    %     D1 = Img_data(S_event(i).PixelIdxList);
    %     mask_D =ismember(D,D1);
    %     D(~mask_D)=nan;
    
    %     Ds = csaps({1:1:size(D,1),1:1:size(D,2)},D,0.1,...
    %         {1:1:size(D,1),1:1:size(D,2)});
    %
    
    
    
    try
    
        % get centroid of event
        centreOfEvent = S_event(i).WeightedCentroid;
        
        centreOfEvent = centreOfEvent - [cols_e(1)-1,rows_e(1)-1];
        
        max_D = max(D(:));
        [r,~] = find(D==max_D,1,'first');
        if abs(centreOfEvent(2)-r)*pxSz > 1
            r = round(centreOfEvent(2));
        end
       
        if (r-(n_px_t-1)/2 <= 0) || (r+(n_px_t-1)/2 > size(D,1))
            
            n_px_t = min(2*r-1,2*(size(D,1)-r)-1);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
%         keyboard
%         
%         
%         % % %%%%%%%
%                    figure
%                    imshow(D,[])
%                    
%                    x = (1:1:size(D,2)).*pxSize_t;
%                    y = (1:1:size(D,1)).*pxSz;
%                    [X,Y] = meshgrid(x,y);
%                    
%                    % [A, pos of peak(x,y),sig_x , sig_y 1]
%                    p0 = [max_D centreOfEvent 5 3 1];
%                    
%                    fun = @(p,x,y) p(1).*exp(-((x-p(2)).^2./(2*p(4).^2))) .* exp(-((y-p(3)).^2./(2*p(5).^2))) + p(6);
%                    
%                    p0 = [max_D centreOfEvent(1) 5 10 centreOfEvent(2) 5 10 1];
%                                      
%                    % coefficients: [A, m1,sd1,tau1, m2,sd2,tau2, F0]
%                    fun = @(p,x,y) p(1) .* exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - x./p(4)).*cdf('Normal',x,p(2)+(p(3)^2/p(4)),p(3)) .* ...   
%                                           exp(p(5)/p(4) + p(6)^2/(2*p(7)^2) - y./p(7)).*cdf('Normal',y,p(5)+(p(6)^2/p(7)),p(6)) + p(8);
%                    
%                    % coefficients: [A, m1,sd1,tau1, m2,sd2, F0]
%                    p0 = [max_D centreOfEvent(1) 5 10 centreOfEvent(2) 5 10 1];
%                    fun = @(p,x,y) p(1) .* exp(p(2)/p(4) + p(3)^2/(2*p(4)^2) - x./p(4)).*cdf('Normal',x,p(2)+(p(3)^2/p(4)),p(3)) .* ...   
%                                           exp(-((y-p(5)).^2./(2*p(6).^2))) + p(8);                   
%                                       
%                    % weights
%                    wE = D./max_D;
%                    
%                    fitFunSum = @(p,x,y,D) sum ( sum( (D-fun(p,x,y)).^2 ) );
%         
%                    opt = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
%                             'MaxIter',1000,'MaxFunEvals',3000,'UseParallel',0);
%                    coef = fmincon(@(p)fitFunSum(p,X,Y,D),p0,[],[],[],[],[],[],[],opt);
%                    
%                    FWHM = 2*sqrt(2*log(2))*coef(6)
%                    
%                    figure                                                                   
%                    surf(X,Y,D,'FaceAlpha',0.5,'FaceColor',[0 0 1],'EdgeColor',[0 0 0])                  
%                    hold on                  
%                    surf(X,Y,fun(coef,X,Y),'FaceAlpha',0.5,'EdgeColor','none','FaceColor',[1 0 0])
% %                 
% %                    xx  = 1
% %                    yy = 2
%         % %          hold on
%         %          plot(centreOfEvent(1),centreOfEvent(2),'ro')
%         %          hold on
%         %          plot(c,r,'go')
% %         
%         figure
%         fnplt(spline1)
%         hold on 
%         fnplt(fnder(spline1))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        t = linspace(0,size(D,2)*pxSize_t,size(D,2));
        t_ups = linspace(0,size(D,2)*pxSize_t,1000);
        
        t_prof = mean(D(r-(n_px_t-1)/2:r+(n_px_t-1)/2,:),1);
        %t_prof_s = smooth(t_prof,n_pts/length(t_prof),'loess');
        
        [v_max_t,p_max_t] = max(t_prof);
        
        % use coeficient from previous fitting (t0, baseline) and calculate
        % tauR
                 
        if ~isempty(prevFitCoef)
            
            t0 = prevFitCoef(i,1) - cols_e(1)*pxSize_t; if t<0, t0=1; end   
            t1 = p_max_t*pxSize_t; 
            A = v_max_t;
            bs1 = prevFitCoef(i,4);
            bs2 = mean(t_prof(end-4:end));
            
        else
            t1 = p_max_t*pxSize_t; 
            t0 = t1-10; if t<0, t0=1; end        
            A = v_max_t;
            bs1 = mean(t_prof(1:5));
            bs2 = mean(t_prof(end-4:end));
            
        end
        
        try
            if t0<t1
                
                [~,indTau] = min(abs(t_prof(round(t0/pxSize_t):p_max_t) - 0.66*A));               
                tR = indTau*pxSize_t;
                
            else
                tR = 5;                
            end
            
        catch
            tR = 5;
            
        end
        
        
        
        % t0, t1, tauR1, tauR2, A, tauD1, tauD2, y01, y02
        %x0 = [t0 p_max_t*pxSize_t 5 5 v_max_t 10 20 1 1];
        
        %          t0, tauR1, A, t1, tauD1, tauD2, y01, y02
        %         x0 = [t0 5 v_max_t p_max_t*pxSize_t 10 20 1 1];
        %keyboard
        
        %         [coeff_t,t_prof_ups,t_prof_fit] = fitWholeSparkParams(t,t_prof,x0,t_ups);
        %
        % %
        % % % %
        %              figure
        %               plot(t,t_prof)
        %               hold on
        %               plot(t,t_prof_fit,'r')
        %
        
        % starting points for chosen fit function
        % coefEst = paramEstimEMG(t_prof,t,t0)        
        % 'EMG' : coeff = [t0, F01, A, m, sd, tau]
        % p0 = [t0 bs1 A t1-tR tR 2*tR];
        
        % 'EGH' : coefficients: [F01, A, m, sd, tau]
        % p0 = [bs1 A t1 tR 2*tR];
             
        % '1expR1expD' : coefficients: [t0, F01, tauR, A, t1, tauD, F02]
        % p0 = [t0 bs1 tR A t1 2*tR bs2];
        
        % '2expR2expD' : coeff = [t0, F01, tauR1, tauR2, A, t1, tauD1, tauD2, F02]
         p0 = [t0 bs1 tR 1.5*tR A t1 2*tR 3*tR bs2];
        
        
        out_t_prof = fitOneWholeEvent(p0,t,t_ups,t_prof,'yes','2expR2expD');
        
%         figure
%         plot(t,t_prof,'ok')
%         hold on 
%         plot(out_t_prof.x_ups,out_t_prof.yFit_ups,'r')
%         
        
        
        
        
%         keyboard
        % coefficients: [t0,F01,A,m,sd,tau,F02]
%         options = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
%             'MaxIter',1000,'MaxFunEvals',3000);
%        coeff_t = fmincon(@(p)expModGauss(p,t,t_prof,'yes','sumOfSquaers'),...
%             [t0 bs v_max_t p_max_t*pxSize_t-5 5 10 1],...
%             [],[],[],[],[],[],[],options);
        
        % t_prof_fit = expModGauss(coeff_t,t,[],'yes','evalFunc');
        % t_prof_ups = expModGauss(coeff_t,t_ups,[],'yes','evalFunc');
        
        t_prof_ups = out_t_prof.yFit_ups;
        coeff_t = out_t_prof.coeff;
        
        %          figure('Name','gauss')
        %          plot(t,t_prof,'k')
        %          hold on
        %          plot(t,t_prof_fit.wholeFit,'r')
        
        [~,p_x] = max(t_prof);
        t_prof_test = t_prof;
        while abs(centreOfEvent(1)-p_x)*pxSize_t > 30
            
            t_prof_test(p_x)=0;
            [~,p_x] = max(t_prof_test);
            
        end
        
        if (p_x-(n_px_x-1)/2 <= 0) || (p_x+(n_px_x-1)/2 > size(D,2))
            
            n_px_x = min(2*p_x-1,2*(size(D,2)-p_x)-1);
            
        end
        
        x = linspace(0,size(D,1)*pxSz-pxSz,size(D,1))';
        x_prof = mean(D(:,p_x-(n_px_x-1)/2:p_x+(n_px_x-1)/2),2);
        
        x_ups = linspace(0,length(x_prof)*pxSz-pxSz,1000);
        
        [v_max_x,p_max_x] = max(x_prof);
  
        % fit x profile with gauss function
        % coefficients: [A,m,sd,tau] + coeff baseline
%         options = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
%             'MaxIter',1000,'MaxFunEvals',3000);
%         coeff_x = fmincon(@(p)expModGauss(p,x,x_prof,'no','sumOfSquaers'),...
%             [v_max_x p_max_x*pxSz-2 1 3 1],...
%             [],[],[],[],[0 0 0 0 0],[inf inf inf inf inf],[],options);
       
        
        out_x_prof = fitOneWholeEvent([v_max_x p_max_x*pxSz-2 1 3 1],x,x_ups,x_prof,'no',[]);
                
        %x_prof_fit = expModGauss(coeff_x,x,[],'no','evalFunc');
        x_prof_ups = out_x_prof.yFit_ups;
   
        




% %                
% figure
% findpeaks(x_prof,x,'SortStr','descend','Annotate','extents');
   
% 

%%%%%%%%
% 
% [FitResults,GOF,baseline,coeff,residual,xi,yi,BootResults]=peakfit([x x_prof],0,0,1,1);
% 
% figure
% plot(x,x_prof)
% hold on
% plot(xi,yi,'r')



%%%%%%%%

% 
% %              [pks_x,locs_x,w_x,p_x] = findpeaks(x_prof,x,'SortStr','descend');
%                  figure    
%                  plot(x,x_prof)
%                      hold on
%                      plot(x_ups,x_prof_ups,'r')
%          
        
        % parameters from time profile fitting        
        baseline = coeff_t(2);
        t0OfSpark = coeff_t(1);
     
        %time parameters
        [val_t,pos_t] = max(t_prof_ups);
        
        half_max_t = (val_t - baseline)/2 + baseline;
        
        y_t1 = t_prof_ups(1:pos_t-1);
        y_t2 = cat(1 ,zeros(pos_t-1,1), t_prof_ups(pos_t:(length(t_prof_ups))));
        
        [~,pos1] = min(abs(y_t1 - half_max_t));
        [~,pos2] = min(abs(y_t2 - half_max_t));
        
        half_max_t_1 = t_ups(pos1);
        half_max_t_2 = t_ups(pos2);
        
        FDHM = half_max_t_2 - half_max_t_1;
        TTP = t_ups(pos_t) - t0OfSpark;
        
        % spatial parameters
        [val_x,pos_x] = max(x_prof_ups);
        
        % somethimes coeff from from fit fall close to zero (bound)
        baseline_x = min(x_prof_ups); 
        

        half_max_x = (val_x(1)-baseline_x)/2 + baseline_x;
        
        y_x1 = x_prof_ups(1:pos_x-1);
        y_x2 = cat(1 ,zeros(pos_x-1,1), x_prof_ups(pos_x:(length(x_prof_ups))));
        
        [~,localMinL] = min(y_x1);
        [~,localMinR] = min(y_x2(pos_x:end));
        localMinR = localMinR + pos_x - 1;
        
        y_x1 = cat(1, zeros(localMinL-1,1), y_x1(localMinL:end));
        y_x2 = cat(1, y_x2(1:localMinR), zeros(length(y_x2)-localMinR,1));
        
        [~,pos1_x] = min(abs(y_x1 - half_max_x));
        [~,pos2_x] = min(abs(y_x2 - half_max_x));
        
        half_max_x_1 = x_ups(pos1_x);
        half_max_x_2 = x_ups(pos2_x);
        
        FWHM = half_max_x_2 - half_max_x_1;
        
        if isempty(FDHM)||isempty(TTP)||isempty(FWHM)
            FDHM = nan;
            TTP = nan;
            FWHM = nan;
        end
        
        Ampl = val_t - baseline;
        
        % calculating mass
        sparkMass = Ampl*1.206*FWHM^3;
        
        eventParams(i).amplitude = Ampl;
        %eventParams(i).duration = max(t);
        eventParams(i).FDHM = FDHM;
        eventParams(i).TTP = TTP;
        eventParams(i).FWHM = FWHM;
        eventParams(i).sparkMass = sparkMass;
        eventParams(i).calciumEventROI = S_event;
        
        %%%%%%%%%%%%%%%%%%%%%%
        %% show  individual calcium events
        if get(getappdata(main_fig,'check_showEventsFigs'),'Value')
            scF = 0.6;
            figure('Tag','CaEventParam','Name',sprintf('posOfRoi %.2f um; event #%d',posOfROI*pxSz,i),'Visible','on',...
                'OuterPosition',[scrRes(3)-min(scrRes(3),scrRes(4))*scF scrRes(3)-min(scrRes(3),scrRes(4))*scF min(scrRes(3),scrRes(4))*scF min(scrRes(3),scrRes(4))*scF])
            subplot(2,2,1)
            image(D,'YData',[min(x_ups) max(x_ups)],'XData',[min(t_ups) max(t_ups)],'CDataMapping','scaled')
            line(get(gca,'XLim'),[(r-(n_px_t-1)/2)*pxSz (r-(n_px_t-1)/2)*pxSz],'Parent',gca,'LineWidth',1,'Color','k',...
                'LineStyle','-');
            line(get(gca,'XLim'),[(r+(n_px_t-1)/2)*pxSz (r+(n_px_t-1)/2)*pxSz],'Parent',gca,'LineWidth',1,'Color','k',...
                'LineStyle','-');
            
            line([(p_x-(n_px_x-1)/2)*pxSize_t (p_x-(n_px_x-1)/2)*pxSize_t],get(gca,'YLim'),'Parent',gca,'LineWidth',1,'Color','k',...
                'LineStyle','-');
            line([(p_x+(n_px_x-1)/2)*pxSize_t (p_x+(n_px_x-1)/2)*pxSize_t],get(gca,'YLim'),'Parent',gca,'LineWidth',1,'Color','k',...
                'LineStyle','-');
            title(sprintf('posOfRoi %.2f \\mum; calcium event (#%d) image ',posOfROI*pxSz,i))
            set(gca,'YDir','normal','FontSize',12)
            xlabel('t (ms)')
            ylabel('x (\mum)')
            
            subplot(2,2,2)
            plot(x_prof,x)
            hold on
            plot(x_prof_ups,x_ups,'r','LineWidth',1)
            set(gca,'FontSize',12)
            title(sprintf('posOfRoi %.2f \\mum; calcium event (#%d) x profile',posOfROI*pxSz,i))
            xlabel('fluorescence (F/F_0)')
            ylabel('x (\mum)')
            ylim([min(x_ups) max(x_ups)])
            line([half_max_x half_max_x],[half_max_x_1 half_max_x_2],'Parent',gca,'LineWidth',1,'Color','k')
            
            subplot(2,2,3)
            plot(t,t_prof)
            hold on
            plot(t_ups,t_prof_ups,'r','LineWidth',1)
            set(gca,'FontSize',12)
            title(sprintf('posOfRoi %.2f \\mum; calcium event (#%d) t profile',posOfROI*pxSz,i))
            xlabel('t (ms)')
            ylabel('fluorescence (F/F_0)')
            xlim([min(t_ups) max(t_ups)])
            line([half_max_t_1 half_max_t_2],[half_max_t half_max_t],'Parent',gca,'LineWidth',1,'Color','k')
            line([t_ups(pos_t) t_ups(pos_t)],[baseline val_t],'Parent',gca,'LineWidth',1,'Color','k')
            line([t0OfSpark t_ups(pos_t)],[baseline baseline],'Parent',gca,'LineWidth',1,'Color','k')
            
            subplot(2,2,4)
            set(gca,'Visible','off')
            text(min(get(gca,'XLim')),max(get(gca,'YLim')/2),...
                {sprintf('posOfRoi %.2f \\mum',posOfROI*pxSz)
                sprintf('parameters of calcium event (#%d):',i)
                sprintf('amplitude = %0.4f (F/F_0)',Ampl)
                sprintf('FDHM = %0.4f (ms)',FDHM)
                sprintf('TTP = %0.4f (ms)',TTP)
                sprintf('FWHM = %0.4f (um)',FWHM)
                sprintf('sparkMass = %0.4f (F/F_0*um^3)',sparkMass)},...
                'Parent',gca,'FontUnits','normalized','FontSize',0.08)
            
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clearvars D t x t_prof x_prof t_ups x_ups t_prof_ups x_prof_ups y_t1 y_t2 D1 y_x1 y_x2
        
    catch ME
        
        eventParams(i).amplitude = nan;
        %eventParams(i).duration = nan;
        eventParams(i).FDHM = nan;
        eventParams(i).TTP = nan;
        eventParams(i).FWHM = nan;
        eventParams(i).sparkMass = nan;
        eventParams(i).calciumEventROI = nan;
        
        
        rethrow(ME)
        
    end
    
    clearvars rows_e cols_e
    
end


end

