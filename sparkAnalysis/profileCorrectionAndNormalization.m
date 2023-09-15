function p_out = profileCorrectionAndNormalization(t,p,main_fig,signalType,typeOfOperation)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% t = time
% p = profile
% signalType = type of signal, 'cyt'(positive peaks) or 'SR' (negative peaks)
% typeOfOperation = 'norm' , 'baseline' or 'both'  


switch typeOfOperation
    
    case 'baseline'
        % baseline correction 
        p_out = bsCorrection(t,p,signalType,main_fig);
                
    case 'norm'
        % normalization of signal 
        p_out = normalization(t,p,signalType);
               
    case 'both'
        % first baseline correction then normalization
        pC = bsCorrection(t,p,signalType,main_fig);
        p_out = normalization(t,pC,signalType);
        
end

end

%%%%% function for baseline correction of signal %%%%%
function pC = bsCorrection(t,p,signalType,main_fig)
% out - pC - signal with corrected baseline
switch signalType
    
    case 'cyt'
        
        %         %%% detred %%%
        %         % calculate break points for detrending, intervals 1000 ms
        %         N_intervals = round(t(end)/1000);
        %         bp = (1:1:N_intervals-1).*1000;
        %         % average of signal in range mean+-2*SD
        %         globalAvrgCYT = mean( p( (p < mean(p)+std(p)) & (p > mean(p)-std(p)) ) );
        %         % detrend, linear
        %         p = detrend(p,'linear',bp) + globalAvrgCYT;
        %
        %         %%% fit %%%
        %         ft = fittype('poly5');
        %         opts = fitoptions( 'Method','LinearLeastSquares',...
        %                             'Normalize','on','Robust','Bisquare');
        %
        %         % fit signal with subtracted trend with high order polynom
        %         % exclude events > 75th percentile
        %         perc = prctile(p,[50 75]);
        %         m_exc = p > perc(2);
        %         opts.Exclude = m_exc;
        %
        %         % fit polynom to signal
        %         [f, ~,~] = fit(t(:),p(:),ft,opts );
        %         pC = feval(f,t);
        pC = p;
        
    case 'SR'
        
        % remove effect of photobleaching
        %%% fit signal with double exp %%%        
        pxSz_t = getappdata(main_fig,'pxSize_t');
                
        ft = fittype('exp2');
        opts = fitoptions('Method','NonLinearLeastSquares',...
            'Normalize','on','Robust','Bisquare');
        [f,~,~] = fit(t(:),p(:),ft,opts);
                
        % remove trend 
        p_dtr = p(:)-feval(f,t(:));
        
        % smoothed signal, moving average, window 50 ms 
        span = ceil(100/pxSz_t)/numel(p);        
        p_dtr_s = smooth(p_dtr(:),span,'rloess');
        
        % calculate median absolute deviation of smoothed signal         
        mad_tr = median(abs(p_dtr_s - median(p_dtr_s)));
        
        % constant, everything less than k*mad_tr considered as blink
        k = 1;
        tresh = median(p_dtr_s)-k*mad_tr;
                       
        % refit signal by excluding signal lower than negative k*mad_tr
        opts.Exclude = p_dtr_s < tresh;
        opts.Lower = [-inf -inf -inf -inf];
        opts.Upper = [inf inf inf inf];
        
        %opts.Weights = dataW;
%           keyboard            
        [f, ~,~] = fit(t(:),p(:),ft,opts);
          
% keyboard
%         figure('Tag','test')
%         plot(t,p)
%         hold on 
%         plot(t,p_dtr,'b')
%         plot(t,feval(f,t),'r')
         
        

         % iterative fitting of baseline, with penalty to data lower than current fit 
%         % old sse     
%         sseO = out.sse;       
%         stopCond = true;
%         % stop when difference between sse is smaller than 5 %
%         while stopCond
%                                   
%             % refit and set weights to points
% %             pridat este pocitanie medianu alebo priemeru a na zaklade toho urcovat vahu bodov 
% %             alebo smoothing a podla toho 
%            
%             dataW = p(:)-feval(f,t);
%             dataW(dataW>0) = 1 - (p(dataW>0)'-(feval(f,t(dataW>0))))./max(p(dataW>0)'-(feval(f,t(dataW>0))));
%             dataW(dataW<0) = p(dataW<0)'./(feval(f,t(dataW<0)).*3);
%             opts.Weights = dataW;
%             
%             [f, out,~] = fit(t(:),p(:),ft,opts);
%                         
%             % new sse
%             sseN = out.sse;
%             
%             % check if diff is smaller than 5% to stop
%             stopCond = ((sseO-sseN)/sseO)*100 > 5;          
%             % update old sse    
%             sseO = sseN;
%             
%         end
%         
%         
%         figure
%         plot(t,p)
%         hold on             
%         plot(t,feval(f,t),'r')
%         hold on 
%         plot(t,p(:)./feval(f,t),'g')
        
        % divide signal wit fit
        pC = p(:)./feval(f,t);
        
end


end

%%%%% function for normalization of signal (self-ratio of signal) %%%%%
function pN = normalization(t,p,signalType)
% out - pC - normalized signal
switch signalType
    % F0 caclulated from whole signal, except of events
    case 'cyt'
        
        F0 = mean(p(p < mean(p) + std(p)));
        
    case 'SR'
     
        %F0 = mean(p(p > mean(p) - 0.5*std(p)));
        F0 = 1;
              
end

pN = p./F0;

end
