function [imgE,imgEs,imgEm,maxOfEventPos,maxCrossProfs,t0,bs] = ...
        eventImgAndMaxCrossingProfiles(...
        statEvent,startOfEvent,endOfEvent,...
        pxSzT,pxSzX,n_px_t,n_px_x,img,prevFitCoef)


% try to expand area around spark in t direction, 
% use start and end of spark from profile
rows_e = statEvent.SubarrayIdx{1};
cols_e = statEvent.SubarrayIdx{2};
% image of event
imgE = img(rows_e,cols_e);

% mask of image
imgMask = false(size(img));
% mask of event in rectangle image of event
imgMask(rows_e,cols_e) = statEvent.Image;

%     % smoothed spark data with spline
%     D = csaps({linspace(1,size(D,1),size(D,1)),linspace(1,size(D,2),size(D,2))},...
%         D,0.25,...
%         {linspace(1,size(D,1),size(D,1)),linspace(1,size(D,2),size(D,2))});

%     figure
%     ax1=subplot(2,1,1)
%     imagesc(D)
%     ax2=subplot(2,1,2)
%     imagesc(Ds)
%     linkprop([ax1,ax2],{'XLim','YLim'})

%     % mask of pixels higher than 90. percentile
%     m_spD = D>prctile(D(:),90);
%     CC_event = bwconncomp(m_spD,8);

%     %calculate properties of all CC (regions) in spark area, to get centre
%     %of spark
%     statOfSubRegions = regionprops(CC_event,D, 'WeightedCentroid','Area','SubarrayIdx');
%
%     % find biggest region
%     [~,p] = max([statOfSubRegions.Area]);

% try to expand spark/event area, but not to have overlaping events in
% time axis
% get position of centre of event
centr = statEvent.WeightedCentroid;
r_m = round(centr(2) - min(rows_e));
c_m = round(centr(1) - min(cols_e));

% dimension for calulation of profiles (t and x) crossing peak
% decrease width of stripe for calculation of time profile
if (r_m-(n_px_t-1)/2 <= 0) || (r_m+(n_px_t-1)/2 > size(imgE,1))
    n_px_t = min(2*r_m-1, abs(2*(size(imgE,1)-r_m)-1));
end

% decrease width of stripe for calculation of spatial profile
if (c_m-(n_px_x-1)/2 <= 0) || (c_m+(n_px_x-1)/2 > size(imgE,2))
    n_px_x = min(2*c_m-1, abs(2*(size(imgE,2)-c_m)-1));
end

% check if isempty startOfEvent, try to estimate it
if isempty(startOfEvent)
    % position of middle of spark in whole image
    r_m_whImg = rows_e(1)-1+r_m;
    % time axis of image
    t_whImg = (1:1:size(img,2)).*pxSzT - pxSzT;
    
    % try to expand area around spark in t direction, fit the rising part
    % of spark, same as in spark recovery analysis
    t_spark_prof_whImg = mean( img(r_m_whImg-(n_px_t-1)/2:r_m_whImg+(n_px_t-1)/2,:) , 1 );
    t_spark_prof_whImg = t_spark_prof_whImg(:);
    
    %fit only rise of spark fun(t0,u,tR,A,y0)
    locs = t_whImg(cols_e(1)-1+c_m);
    pks = t_spark_prof_whImg(cols_e(1)-1+c_m);
    
    try
        [~,~,prevFitCoef,~,startOfEvent,endOfEvent] = fitSparkRise(pxSzT,t_whImg,...
            t_spark_prof_whImg,pks,locs,[],[],1e-9,1000,50,75,[],[]);
        
    catch
        prevFitCoef = [];
        startOfEvent = [];
        endOfEvent = [];
    end
    
end

if ~isempty(startOfEvent)
    s = startOfEvent;
    e = endOfEvent;
    if s>cols_e(1), s = cols_e(1); end
    if e<cols_e(end), e = cols_e(end); end
    % adjusted position of center
    c_m = c_m + (cols_e(1)-s);
    cols_e = s:e;
end

% adjusted event area data
imgE = img(rows_e,cols_e);
imgEm = imgMask(rows_e,cols_e);

% % smooth little bit, to have only one max
% imgEs = csaps({linspace(1,size(imgE,1),size(imgE,1)),linspace(1,size(imgE,2),size(imgE,2))},...
%     imgE,0.5,...
%     {linspace(1,size(imgE,1),size(imgE,1)),linspace(1,size(imgE,2),size(imgE,2))});
% 
% % event data only, else NaN
% %D(~D_m) = nan;
% imgEs(~imgE_m) = nan;

%r_m = round(centr(2)) - rows_e(1)+1;
%[r_m,c_m] = find(imgEs==max(imgEs(:)));

% % remove areas where max crossing profile is nan
% pT = mean( imgEs( r_m-(n_px_t-1)/2:r_m+(n_px_t-1)/2, : ), 1 );
% pX = mean( imgEs( :, c_m-(n_px_x-1)/2:c_m+(n_px_x-1)/2 ), 2 );
% 
% pT_m = isnan(pT);
% pX_m = isnan(pX);
% 
% % check that there are non NaNs somwhere in middle of profile
% ind_beforeMax_t = find(pT_m(1:c_m),1,'last');
% ind_afterMax_t = find(pT_m(c_m:end),1,'first') + c_m - 1;
% 
% ind_beforeMax_x = find(pX_m(1:r_m),1,'last');
% ind_afterMax_x = find(pX_m(r_m:end),1,'first') + r_m - 1;
% 
% if ~isempty(ind_beforeMax_t)
%     pT_m(1:ind_beforeMax_t) = true;
% end
% 
% if ~isempty(ind_afterMax_t)
%     pT_m(ind_afterMax_t:end) = true;
% end
% 
% if ~isempty(ind_beforeMax_x)
%     pX_m(1:ind_beforeMax_x) = true;
% end
% 
% if ~isempty(ind_afterMax_x)
%     pX_m(ind_afterMax_x:end) = true;
% end
% 
% % adjusted event area data
% imgE = img(rows_e(~pX_m),cols_e(~pT_m));
% imgE_m = imgMask(rows_e(~pX_m),cols_e(~pT_m));


%% redo smoothing and width of profiles estimation 
% smooth little bit, to have only one max
try
    imgEs = csaps(...
        {linspace(1,size(imgE,1),size(imgE,1)),linspace(1,size(imgE,2),size(imgE,2))},...
        imgE,0.5,...
        {linspace(1,size(imgE,1),size(imgE,1)),linspace(1,size(imgE,2),size(imgE,2))});
    
    % event data only, else NaN
    %D(~D_m) = nan;
    imgE_s = imgEs;
    imgE_s(~imgEm) = nan;
    
    % get position of centre of event
    %centr = statEvents(i).WeightedCentroid;
    %r_m = round(centr(2)) - rows_e(1)+1;
    [r_m,c_m] = find(imgE_s==max(imgE_s(:)));
    
%     keyboard
%     figure
%     imagesc(imgEs)

    % dimension for calulation of profiles (t and x) crossing peak
    % decrease width of stripe for calculation of time profile
    if (r_m-(n_px_t-1)/2 <= 0) || (r_m+(n_px_t-1)/2 > size(imgE,1))
        n_px_t = min(2*r_m-1, abs(2*(size(imgE,1)-r_m)-1));
    end
    
    % decrease width of stripe for calculation of spatial profile
    if (c_m-(n_px_x-1)/2 <= 0) || (c_m+(n_px_x-1)/2 > size(imgE,2))
        n_px_x = min(2*c_m-1, abs(2*(size(imgE,2)-c_m)-1));
    end
    
    [~,c_m] = max( mean( imgE( r_m-(n_px_t-1)/2:r_m+(n_px_t-1)/2, : ), 1 ) );
    %c_m = c_m + statOfSubRegions(p).SubarrayIdx{2}(1) - 1;
    
    
    % create time and spatial axis and profiles
    t = (1:1:size(imgE,2)).*pxSzT - pxSzT; % x axes = time
    x = (1:1:size(imgE,1)).*pxSzX - pxSzX; % y axes = spatial
    t = t(:);
    x = x(:);
    
    % upscale 10x
    uspFactor = 10;
    t_ups = linspace(t(1),t(end),numel(t)*uspFactor); % x axes = time
    x_ups = linspace(x(1),x(end),numel(x)*uspFactor); % y axes = spatial
    t_ups = t_ups(:);
    x_ups = x_ups(:);
    
    % x and t max crossing profiles
    t_event_prof = mean( imgE(r_m-(n_px_t-1)/2:r_m+(n_px_t-1)/2,:) , 1 );
    t_event_prof = t_event_prof(:);
    x_event_prof = mean( imgE(:,c_m-(n_px_x-1)/2:c_m+(n_px_x-1)/2) , 2 );
    x_event_prof = x_event_prof(:);
    
    % x and t max crossing profiles from smoothed event
    t_event_prof_S = mean( imgEs(r_m-(n_px_t-1)/2:r_m+(n_px_t-1)/2,:) , 1 );
    t_event_prof_S = t_event_prof_S(:);
    x_event_prof_S = mean( imgEs(:,c_m-(n_px_x-1)/2:c_m+(n_px_x-1)/2) , 2 );
    x_event_prof_S = x_event_prof_S(:);
    
    % % look 10 ms around (in points) column with detected event maximum
    % eps = ceil(10/pxSzT);
    % startPoint = c_m-eps;
    % if startPoint<1, startPoint=1; end
    % endPoint = c_m+eps;
    % if endPoint>numel(t_spark_prof), endPoint=numel(t_spark_prof); end
    %
    % [v_max_t,p_max_t,~,~] = findpeaks( ...
    %     t_spark_prof(startPoint:endPoint),...
    %     'SortStr','descend');
    % p_max_t = p_max_t + startPoint - 1;
    %
    % peakData.pos = p_max_t;
    % peakData.val = v_max_t;
    
    [v_max_t,p_max_t] = max(t_event_prof);
    peakData_tProf.pos = p_max_t;
    peakData_tProf.val = v_max_t;
    
    % use coeficient from previous fitting (t0, tauR, A, baseline) to
    % estimate bs, t0
    if ~isempty(prevFitCoef)
        try
            t0 = prevFitCoef(1) - cols_e(1)*pxSzT;
            bs = mean(t_event_prof(1:ceil(t0/pxSzT)));
        catch
            t0 = p_max_t*pxSzT - 5; % in ms, 5 ms before max
            if t0<1, t0 = 1; end
            bs = mean(t_event_prof(1:ceil(t0/pxSzT)));
        end
    else
        t0 = p_max_t*pxSzT - 5; % in ms, 5 ms before max
        if t0<1, t0 = 1; end
        bs = mean(t_event_prof(1:ceil(t0/pxSzT)));
    end
    
catch
    %% too few points, failed
    imgE = ones([10,10]);
    imgEm = true(size(imgE));
    imgEs = ones([10,10]);
    r_m = 1;
    c_m = 1;
    t = ones([1,10]);
    t_ups = ones([1,10]);
    t_event_prof = ones([1,10]);
    t_event_prof_S = ones([1,10]);
    peakData_tProf = [1,1];
    x = ones([1,10]);
    x_ups = ones([1,10]);
    x_event_prof = ones([1,10]);
    x_event_prof_S = ones([1,10]);
    t0 = 1;
    bs = 1;
    
    
end
 

%% output
maxOfEventPos = [r_m,c_m];

maxCrossProfs.t = t;
maxCrossProfs.t_ups = t_ups;
maxCrossProfs.t_event_prof = t_event_prof;
maxCrossProfs.t_event_prof_S = t_event_prof_S;
maxCrossProfs.peakData_tProf = peakData_tProf;
maxCrossProfs.x = x;
maxCrossProfs.x_ups = x_ups;
maxCrossProfs.x_event_prof = x_event_prof;
maxCrossProfs.x_event_prof_S = x_event_prof_S;


end

