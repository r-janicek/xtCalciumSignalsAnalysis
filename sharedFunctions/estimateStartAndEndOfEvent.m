function [pos_s, pos_e] = estimateStartAndEndOfEvent( ...
    prof, peakPosPx, options)
%{
estimate start and end of event with maximum defined baseline
params:
prof = fluorescence profile
peakPosPx = position of peak
options:
    maxDurOfBaseline = maximum duration of baseline
    equalBaselineDur = symetric maximum baseline
    evntsMask = mask of event in profile
    smoothSpan = number of points to consider 
    evntAcceptCrit = in percentage
%}
arguments         
    prof (:,1) {mustBeNumeric}
    peakPosPx (1,1) {mustBeInteger, mustBePositive}
    options.maxDurOfBaseline (1,1) {mustBeInteger, mustBePositive} = 10
    options.equalBaselineDur (1,1) logical = true
    options.evntsMask (:,1) logical = true(size(x))
    options.smoothSpan (1,1) {mustBeInteger, mustBePositive} = 10
    options.evntAcceptCrit (1,1) {mustBeInteger, mustBePositive} = 75
    options.baselineDurMult (1,1) {mustBeNumeric, mustBePositive} = 2
end

% find start and end of event with baseline
% mask of fitted event
if ~all(options.evntsMask, 'all')
    try
        % just in case there are also other events masked
        m_event = false(size(prof));
        evnt_m_s = find(~options.evntsMask(1:peakPosPx), 1, 'last');
        if isempty(evnt_m_s), evnt_m_s = 1; end
        evnt_m_e_afterPeak = find(~options.evntsMask(peakPosPx:end), 1, 'first')-1;
        if isempty(evnt_m_e_afterPeak), evnt_m_e_afterPeak = 0; end
        evnt_m_e = peakPosPx + evnt_m_e_afterPeak;
        m_event(evnt_m_s:evnt_m_e) = true;
    catch
        m_event = options.evntsMask;
    end
else
    m_event = options.evntsMask;
end

% take event peak position and its surrounding with maximum baseline 
m_eventWithBsl = false(size(prof));
m_eventWithBsl( ...
    max(1, peakPosPx-options.maxDurOfBaseline) : ...
    min(numel(prof), peakPosPx+ceil(options.baselineDurMult*options.maxDurOfBaseline))) = true;
% smooth event with loess to estimate 
% remove baseline, only events bigger than specified percentile stay
prof_s = nan(size(prof));
prof_s(m_eventWithBsl) = smooth(prof(m_eventWithBsl), ...
    min(0.99, options.smoothSpan/numel(prof(m_eventWithBsl))), ...
    'loess');
% calculate threshold from smoothed event
percentl = prctile(prof_s(m_eventWithBsl & prof_s>0), ...
    [25 50 options.evntAcceptCrit]);
if max(percentl) > max(prof_s(m_event))
    percentl = prctile(prof_s(m_event), ...
        [25 50 options.evntAcceptCrit]);
end

%iqr = percentl(3)- percentl(1);
bsl_thrsh = percentl(3);
prof_s(isnan(prof_s)) = bsl_thrsh;
% treshold profile
prof_s(prof_s < bsl_thrsh) = bsl_thrsh;
% get all posible peaks in smoothed profile of events
try
    [valPeaks_s, locPeaks_s] = findpeaks(prof_s(m_event));
    if isempty(locPeaks_s)
        [valPeaks_s, locPeaks_s] = max(prof_s(m_event));
    end
catch
    [valPeaks_s, locPeaks_s] = max(prof_s(m_event));
end
locPeaks_s = locPeaks_s + find(m_event, 1, 'first') - 1;
% get the one closest to peak of event
[~, idx_p] = min(abs(locPeaks_s-peakPosPx));
% change values of currently fitted event to peak value
% (to be sure I have correct estimate of baseline.
% It might happen that detected spark has multiple peaks in decaying phase,
% this will remove them and still keep nice trace to calculate gradient)
% check
m_event_thrs = prof_s > bsl_thrsh;
if numel(locPeaks_s) > 1
    prof_s_afterPeak = valPeaks_s(idx_p)*m_event_thrs(peakPosPx:end);
    prof_s_afterPeak(prof_s_afterPeak==0) = bsl_thrsh;
    prof_s(peakPosPx:end) = prof_s_afterPeak;
end

% calculate start and end of event using derivation
% calculated on thresholded event profile
% find start of event
prof_s_beforePeak = prof_s(1:locPeaks_s(idx_p));
% flip profile to calulate gradient from left to right
prof_s_beforePeak = flipud(prof_s_beforePeak(:));
pos_s = max( [...
    locPeaks_s(idx_p)-find(gradient(prof_s_beforePeak)>0, 1, 'first')+1, ...
    locPeaks_s(idx_p)-options.maxDurOfBaseline ] );
if isempty(pos_s) || isnan(pos_s) || pos_s < 1
    pos_s = find(m_eventWithBsl, 1, 'first');
    if isempty(pos_s)
        pos_s = 1;
    end
end
pos_s = round(pos_s);

% find end of event
prof_s_afterPeak = prof_s(locPeaks_s(idx_p):end);
if options.equalBaselineDur
    % symetric baseline
    pos_e = min( [...
        locPeaks_s(idx_p)+find(gradient(prof_s_afterPeak)>0, 1, 'first')-1, ...
        locPeaks_s(idx_p)+options.maxDurOfBaseline] );
else
    pos_e = min( [...
        locPeaks_s(idx_p)+find(gradient(prof_s_afterPeak)>0, 1, 'first')-1, ...
        locPeaks_s(idx_p)+ceil(options.baselineDurMult*options.maxDurOfBaseline)] );
end
if isempty(pos_e) || isnan(pos_e) || pos_e>numel(prof)
    pos_e = find(m_eventWithBsl, 1, 'last');
end
if pos_e<=pos_s, pos_e = pos_s+1; end
pos_e = round(pos_e);

end



