function [imgE, imgEs, imgEm, maxOfEventPos, ...
    maxCrossProfs, t0, bs, eventROIstart_t] = ...
        eventImgAndMaxCrossingProfiles(...
        statEvent, startOfEvent, endOfEvent,...
        pxSzT, pxSzX, n_px_t, n_px_x, img, prevFitCoef, ...
        bsDetSensitivity, smoothSpan)

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
% try to update WeightedCentroid position
% there should be only one event of interest in imgE mask
% do threshold 90th percentile, create mask and do AND operation with
% previous event mask
%imgE_s = imgaussfilt(imgE, 1);
imgE_trsh_m = imgE >= prctile(imgE(statEvent.Image), 90);
imgE_trsh_m = imgE_trsh_m & statEvent.Image;

% imgE_trsh_m = statEvent.Image;

% calculate properties of all CC (regions) in event area, to get center
% of event
CC_SubRegions = bwconncomp(imgE_trsh_m,8);
statOfSubRegions = regionprops(CC_SubRegions, imgE, ...
    'WeightedCentroid', 'Centroid', 'Area', 'SubarrayIdx', ...
    'MeanIntensity', 'MaxIntensity', 'PixelIdxList', ...
    'Solidity', 'EulerNumber');
% find biggest region
% statOfSubRegions = struct2table(statOfSubRegions, 'AsArray',true);
% statOfSubRegions = sortrows(statOfSubRegions, 'Area', 'descend');

[~,p1] = max([statOfSubRegions.Area]);
[~,p2] = max([statOfSubRegions.MeanIntensity]);
[~,p3] = max([statOfSubRegions.MaxIntensity].*[statOfSubRegions.Area]);
p = mode([p1,p2,p3]);

% leave only selected region in mask
imgE_trsh_m = false(size(imgE_trsh_m));
imgE_trsh_m(statOfSubRegions(p).PixelIdxList) = true;

% figure
% tiledlayout(2,1)
% nexttile
% imagesc(imgE)
% nexttile
% imagesc(imgE_trsh_m)

% get position of center of event
r_m = round(statOfSubRegions(p).WeightedCentroid(2));
% c_m = round(statOfSubRegions(p).WeightedCentroid(1));
if r_m > size(imgE, 1)
    r_m = round(statOfSubRegions(p).Centroid(2));
end

% % get position of center of event
% centr = statEvent.WeightedCentroid;
% r_m = round(centr(2) - min(rows_e));
% c_m = round(centr(1) - min(cols_e));

% dimension for calculation of profiles (t and x) crossing peak
% decrease width of stripe for calculation of time profile
if (r_m-(n_px_t-1)/2 <= 0) || (r_m+(n_px_t-1)/2 > size(imgE,1))
    n_px_t = min(2*r_m-1, abs(2*(size(imgE,1)-r_m)-1));
end
% find column as max in profile of event
prof_t_imgE = mean(imgE(r_m-(n_px_t-1)/2:r_m+(n_px_t-1)/2,:), 1);
prof_t_imgE_m = any(imgE_trsh_m(r_m-(n_px_t-1)/2:r_m+(n_px_t-1)/2,:), 1);

% set values outside of detected region around event maximum as -INF
prof_t_imgE(~prof_t_imgE_m) = -inf; 
[~,c_m] = max(prof_t_imgE);

% decrease width of stripe for calculation of spatial profile
if (c_m-(n_px_x-1)/2 <= 0) || (c_m+(n_px_x-1)/2 > size(imgE,2))
    n_px_x = min(2*c_m-1, abs(2*(size(imgE,2)-c_m)-1));
end

% position of middle of spark in whole image
r_m_whImg = rows_e(1)-1+r_m;
% time axis of image
t_whImg = (1:1:size(img,2)).*pxSzT - pxSzT;
% try to expand area around spark in t direction, fit the rising part
% of spark, same as in spark recovery analysis
t_spark_prof_whImg = mean( ...
    img(r_m_whImg-(n_px_t-1)/2:r_m_whImg+(n_px_t-1)/2,:) , 1 );
t_spark_prof_whImg = t_spark_prof_whImg(:);
t_whImg_evnt_m = any( ...
    imgMask(r_m_whImg-(n_px_t-1)/2:r_m_whImg+(n_px_t-1)/2,:) , 1 );
% check if isempty startOfEvent, try to estimate it
if isempty(startOfEvent)
    % fit only rise of spark fun(t0,tR,A,y0)
    locs = t_whImg(cols_e(1)-1+c_m);
    pks = t_spark_prof_whImg(cols_e(1)-1+c_m);
    if isempty(prevFitCoef)
        try    
            evntRiseFit = fitEventRise(t_whImg, t_spark_prof_whImg, ...
                "local", peaks_vals=pks, peaks_locs=locs, ...
                smooth_span=smoothSpan, bs_crit=bsDetSensitivity, ...
                evntsMask=t_whImg_evnt_m);
            prevFitCoef = evntRiseFit.coef;
            startOfEvent = evntRiseFit.startOfEvent;
            endOfEvent = evntRiseFit.endOfEvent;
        catch
            prevFitCoef = [];
            startOfEvent = cols_e(1);
            endOfEvent = cols_e(end);
        end
    end
end
if ~isempty(startOfEvent)
    s = startOfEvent;
    e = endOfEvent;
    if s>cols_e(1), s = cols_e(1); end
    if e<cols_e(end), e = cols_e(end); end
    dc = abs(s - cols_e(1));
    cols_e = s:e;
end
% mask of t profile of event, take in to account also fit with exp. rise
% function
t_event_prof_m = t_whImg_evnt_m(cols_e);%false(size(cols_e));
% t_event_prof_m(startOfEvent-cols_e(1)+1:endOfEvent-cols_e(1)+1) = true;

% expand spark area in x dimension
[~, c_m_t_bs] = max(t_spark_prof_whImg(cols_e));
c_m_whImg = cols_e(1)-1+c_m_t_bs;
x_spark_prof_whImg = ...
    mean( img(:, max(1,c_m_whImg-(n_px_x-1)/2):min(size(img,2),c_m_whImg+(n_px_x-1)/2)), 2);
x_spark_prof_whImg_m = false(size(x_spark_prof_whImg));
x_spark_prof_whImg_m(rows_e) = true;
% position of peak in x profile
if all(x_spark_prof_whImg_m)
    [~, peakXPosPx] = max(x_spark_prof_whImg); 
else
    x_spark_prof_whImg_evnt = x_spark_prof_whImg;
    x_spark_prof_whImg_evnt(~x_spark_prof_whImg_m) = -Inf;
    [~, peakXPosPx] = max(x_spark_prof_whImg_evnt);
end
% try to expand in x dimension
try
    % use maximum 540 um baseline
    [startOfEvntX, endOfEvntX] = estimateStartAndEndOfEvent( ...
        x_spark_prof_whImg, peakXPosPx, ...
        maxDurOfBaseline=ceil(5/pxSzX), ...
        evntsMask=x_spark_prof_whImg_m, ...
        equalBaselineDur=true, ...
        smoothSpan=round(3/pxSzX), ...
        evntAcceptCrit=bsDetSensitivity);
catch
    startOfEvntX = rows_e(1);
    endOfEvntX = rows_e(end);
end
if ~isempty(startOfEvntX)
    if startOfEvntX>rows_e(1), startOfEvntX = rows_e(1); end
    if endOfEvntX<rows_e(end), endOfEvntX = rows_e(end); end
    dc_x = abs(startOfEvntX - rows_e(1));
    rows_e = startOfEvntX:endOfEvntX;
end

% adjusted event area data
imgE = img(rows_e,cols_e); 


%imgEm = imgMask(rows_e,cols_e);
imgE_subRegion_m = false(size(imgMask));
imgE_subRegion_m(statOfSubRegions(p).SubarrayIdx{1} + rows_e(1) - 1 + dc_x, ...
    statOfSubRegions(p).SubarrayIdx{2} + cols_e(1) - 1 + dc) = true;
imgEm = imgE_subRegion_m(rows_e, cols_e);

eventROIstart_t = min(cols_e(:))*pxSzT;
if isempty(eventROIstart_t), eventROIstart_t = 0; end

%% redo smoothing and width of profiles estimation 
% smooth event image little bit, to have only one max
try 
    imgEs = csaps(...
        {linspace(1,size(imgE,1),size(imgE,1)),linspace(1,size(imgE,2),size(imgE,2))},...
        imgE,0.5,...
        {linspace(1,size(imgE,1),size(imgE,1)),linspace(1,size(imgE,2),size(imgE,2))});
    % event data only, else NaN
    imgE_s = imgEs;
    imgE_s(~imgEm) = nan;
    % get position of centre of expanded event
    [r_m, ~] = find(imgE_s==max(imgE_s,[],'all','omitmissing') );
    if (r_m-(n_px_t-1)/2 <= 0) || (r_m+(n_px_t-1)/2 > size(imgE_s,1))
        n_px_t = min(2*r_m-1, abs(2*(size(imgE_s,1)-r_m)-1));
    end
    prof_t_imgE_s = mean(imgE_s(r_m-(n_px_t-1)/2:r_m+(n_px_t-1)/2,:), ...
        1, "omitmissing");
    %[r_m,c_m] = find(imgE_s==max(imgE_s(:)));
    [~, c_m] = max(prof_t_imgE_s);
    
    % dimension for calulation of profiles (t and x) crossing peak
    % decrease width of stripe for calculation of time profile
    if (r_m-(n_px_t-1)/2 <= 0) || (r_m+(n_px_t-1)/2 > size(imgE,1))
        n_px_t = min(2*r_m-1, abs(2*(size(imgE,1)-r_m)-1));
    end
    
    % decrease width of stripe for calculation of spatial profile
    if (c_m-(n_px_x-1)/2 <= 0) || (c_m+(n_px_x-1)/2 > size(imgE,2))
        n_px_x = min(2*c_m-1, abs(2*(size(imgE,2)-c_m)-1));
    end
    
    % find maximum of time profile (smoothed)
    [~,c_m] = max( ...
        mean(imgE_s(r_m-(n_px_t-1)/2:r_m+(n_px_t-1)/2, :), 1, "omitnan" ) );
    %c_m = c_m + statOfSubRegions(p).SubarrayIdx{2}(1) - 1;

    % look at gradient descent from maximum point to increase mask of event
    profDesc = mean(imgEs(r_m-(n_px_t-1)/2:r_m+(n_px_t-1)/2,1:c_m), 1, ...
        "omitnan" );
    
    sEventEstGrad = c_m - find(gradient(fliplr(profDesc))>0,1,'first') + 1;
    if isempty(sEventEstGrad)
        sEventEstGrad=1;
    elseif sEventEstGrad<1 
        sEventEstGrad=1;
    end
         
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
    
    % peak data of time profile, use mask
    t_event_prof_EventMask = t_event_prof;
    %t_event_prof_EventMask(~imgEm(r_m,:)) = nan;
    t_event_prof_EventMask(~t_event_prof_m) = nan;
    
    [v_max_t,p_max_t] = max(t_event_prof_EventMask);
    peakData_tProf.pos = p_max_t;
    peakData_tProf.val = v_max_t;
    
    % use coeficient from previous fitting (t0, tauR, A, baseline) to
    % estimate bs, t0
    if ~isempty(prevFitCoef)
        try
            t0 = prevFitCoef(1) - cols_e(1)*pxSzT;
            if t0<1, t0 = 1; end
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
    
    % check if t0 is lower than position of maximum, also bs is smaller
    % than maximum value
    if t0>t(p_max_t)
        t0 = t(sEventEstGrad);
    end
    if bs >= v_max_t
        bs = t_event_prof(sEventEstGrad);
    end
  
catch
    %keyboard
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
    t_event_prof_m = ones([1,10]);
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
maxCrossProfs.t_event_prof_m = t_event_prof_m;
maxCrossProfs.peakData_tProf = peakData_tProf;
maxCrossProfs.x = x;
maxCrossProfs.x_ups = x_ups;
maxCrossProfs.x_event_prof = x_event_prof;
maxCrossProfs.x_event_prof_S = x_event_prof_S;


end

