function params = calcParametersOfEventsFromWholeProfileFit(x_t,fitData)

% calculate parameters [t0 (start of event), tPeak (time pos of peak), 
% Apeak (amplitude of peak), TTP (time to peak), FDHM (full duration in half maximum)] of events from whole profile fit 

keyboard
if isfield(fitData,'t_ups')
    % upscaled time axes
    allEventsFits = cell2mat(fitData.t_ups.individualEventsFits(2:end,:));
    t = fitData.t_ups.t_ups;
else
    allEventsFits = cell2mat(fitData.individualEventsFits(2:end,:));
    t = x_t;  
end

% exp decay function 
funExpFit = @(p,x,y) p(1)*exp(-(x-p(2))/p(3)) - y;
funExpEval = @(p,x) p(1)*exp(-(x-p(2))/p(3));

eventsParams = zeros(size(allEventsFits,2),9);
for i = 1:size(allEventsFits,2)
    
    % single event
    eventFit = allEventsFits(:,i);
    
    % calculate time params of single event
    % position and value of peak
    [val_t,pos_t] = max(eventFit);
    
    % value of half max, in this case baseline is zero
    half_max = val_t/2;
    
    % split fit in two parts (peak position) to find positions of half maxima
    y_t1 = eventFit(1:pos_t-1);
    y_t2 = cat(1 ,zeros(pos_t-1,1), eventFit(pos_t:(length(eventFit))));
    
    [~,pos1] = min(abs(y_t1 - half_max));
    [~,pos2] = min(abs(y_t2 - half_max));
    
    half_max_1 = t(pos1);
    half_max_2 = t(pos2);
         
    % find t0, start of spark, 0.05 set manually, not really nice
    %t0 = t(find(eventFit>max(eventFit)/50,1,'first'));
    
    % find t0 (start of spark) as a cross section of 0 and line fit of part of event with amplitudes between 25 and 75 % of max amplitude
    % find position of 25 and 75% of amplitude
    indx_25 = (find(eventFit>max(eventFit)*0.25,1,'first'));
    indx_75 = (find(eventFit>max(eventFit)*0.75,1,'first'));
    
    if indx_25 == indx_75
        
        indx_25 = indx_25 - 1;
        
    end
    
    try
        f_line25_75 = fit(t(indx_25:indx_75),eventFit(indx_25:indx_75),'poly1');
        t_fit = t(find(eventFit>1e-9,1,'first'):pos_t);
        line25_75 = feval(f_line25_75,t_fit);
        t0 = t_fit(find(line25_75<0,1,'last'));
    catch
        t0 = t(find(eventFit>max(eventFit)/50,1,'first'));
    end
    
    
%     figure
%     plot(t,eventFit,'ok')
%     hold on
%     plot(t(indx_25:indx_75),eventFit(indx_25:indx_75),'r')
%     hold on
%     plot(t_fit,line25_75,'b')
    
    %     figure
    %     plot(fitData.wholeFit)
    
    if isempty(t0)
        % find t0, start of spark, 0.05 set manually, not really nice
        t0 = t(find(eventFit>max(eventFit)/50,1,'first'));
    end
           
    % calculate parameters
    FDHM = half_max_2 - half_max_1;
    TTP = t(pos_t) - t0;
    
    % calculate area under curve
    AUC = trapz(t,eventFit);
    
    % calculate 1. derivative max amplitude and position
    dt = mean(diff(t));
    [firstDerVal,firstDerPos] = max(gradient(eventFit)./dt); 
       
    % get tau of decay
    try
 
        ind = strcmp(fitData.coefficientsOfFittedEvents(1,:),'tauD');
        try 
            tauD = fitData.coefficientsOfFittedEvents{i+1,ind};
        catch 
            tauD = [];
        end
        
        ind_r = strcmp(fitData.coefficientsOfFittedEvents(1,:),'tauR');
        try 
            tauR = fitData.coefficientsOfFittedEvents{i+1,ind_r};
        catch 
            tauR = [];
        end
        
        if ~isempty(tauR) && ~isempty(tauD) && tauR > tauD
            tauD = [];
        end
        
        if isempty(tauD)
            % fit with single exponential             
           t_exp = fitData.t;
           y_exp = cell2mat(fitData.individualEventsFits(2:end,i));
           [vM,pM] = max(y_exp);
       
           coefExp = lsqnonlin(@(p)funExpFit(p,t_exp(pM:end),y_exp(pM:end)),[vM,t_exp(pM),20],[0 0 0],[inf inf inf]);
           
           tauD = coefExp(3);
           
%            figure
%            plot(t_exp,y_exp,'ok')
%            hold on
%            plot(t_exp(pM:end),funExpEval(coefExp,t_exp(pM:end)),'r')
           
        end
               
    catch
        tauD = nan;                       
    end
  
    % final matrix of parameters
    try
        eventsParams(i,:) = [t0, t(pos_t), val_t, TTP, FDHM, AUC, tauD, firstDerVal, t(firstDerPos)];
    catch
        eventsParams(i,:) = [nan, nan, nan, nan, nan, nan, nan, nan, nan];
    end
    
    clearvars y_t1 y_t2
    
end

% save data 
names_params = {'t0 (ms)','tPeak (ms)','Apeak (deltaF/F0)','TTP (ms)','FDHM (ms)','AUC (ms*deltaF/F0)','tauD (ms)','firstDerMaxVal (deltaF/F0)','firstDerMaxValPos (ms)'};
params = [names_params;num2cell(eventsParams)];


end

