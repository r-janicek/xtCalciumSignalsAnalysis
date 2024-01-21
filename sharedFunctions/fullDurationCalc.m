function [FD, prcOfAmpl_val, pos_t_1, pos_t_2, pos_1, pos_2] = ...
    fullDurationCalc(t, t_prof, sE, eE, peakVal, peakPos, bs, prcOfAmpl)
% calculate full duration at specified percentage of amplitude of event
% half maximum value
% t = time
% t_prof = profile
% sE, eE = start and end of event
% peakVal, peakPos, bs = values of peak and baseline
% prcOfAmpl = percentage of amplitude of event where to calculate full
% duration
t_indx = (1:1:numel(t));
t_indx = t_indx(:);

% value of desired percentage of amplitude
prcOfAmpl_val = (peakVal - bs)*(prcOfAmpl/100) + bs;

% parts of profiles before and after maximum of event
t_prof_beforePeak = t_prof(1:peakPos-1);
t_prof_beforePeak = flipud(t_prof_beforePeak(:));

% before peak
sE = min(sE, ...
    numel(t_prof_beforePeak) - ...
    find((t_prof_beforePeak-prcOfAmpl_val) < 0, 1, 'first'));
% after peak
t_prof_afterPeak = t_prof(peakPos:end);
p_prcOfAmpl_afterPeak = peakPos + ...
    find((t_prof_afterPeak - prcOfAmpl_val) < 0, 1, 'first');
if (isempty(p_prcOfAmpl_afterPeak) || p_prcOfAmpl_afterPeak>numel(t_prof))
    p_prcOfAmpl_afterPeak = numel(t_prof);
end

eE = max(eE, p_prcOfAmpl_afterPeak);

if sE>1
    y_t1 = [ zeros(sE-1,1); t_prof(sE:peakPos-1) ];
else
    y_t1 = t_prof(1:peakPos-1);
end
y_t2 = [ zeros(peakPos-1,1); t_prof(peakPos:eE) ];

% find position of first half max position
y_t1_r = flipud(y_t1);
d_y_t1 = gradient(y_t1_r);

if all(d_y_t1==0)
    pos_1 = numel(d_y_t1);
else
    indDer1 = find( d_y_t1<=0 & y_t1_r<prcOfAmpl_val, 1, 'first');
    if ~isempty(indDer1), y_t1_r(indDer1:end) = 0; end

    y_t1 = flipud(y_t1_r);
    pos_last_0_y_t1 = find(y_t1>0, 1, 'first')-1;
    % get the first position of desired percentage of amplitude
    [~,pos_1] = min(abs(y_t1(pos_last_0_y_t1:end) - prcOfAmpl_val));
    pos_1 = pos_1 + pos_last_0_y_t1 - 1;
end

% find the second position of desired percentage of amplitude
d_y_t2 = gradient(y_t2);

indDer2 = find( d_y_t2<0 & y_t2<prcOfAmpl_val & t_indx(1:eE)>peakPos, ...
    1, 'first');
if ~isempty(indDer2), y_t2(indDer2:end) = 0; end

pos_last_numNot0_y_t2 = find(y_t2(peakPos:end)>0, 1, 'last') + peakPos - 1;
[~, pos_2] = ...
    min(abs(y_t2(peakPos:pos_last_numNot0_y_t2) - prcOfAmpl_val));
pos_2 = pos_2 + peakPos - 1;

% positions in time
pos_t_1 = t(pos_1);
pos_t_2 = t(pos_2);
% get full duration in half maximum
FD = pos_t_2 - pos_t_1;
if FD < 0
    FD = 0;
elseif isempty(FD)
    FD = 0;
end

end