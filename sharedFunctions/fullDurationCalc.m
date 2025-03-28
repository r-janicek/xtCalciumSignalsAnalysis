function [FD, prcOfAmpl_val, pos_t_1, pos_t_2, pos_1, pos_2] = ...
    fullDurationCalc(t, t_prof, sE, eE, peakVal, peakPos, bs, prcOfAmpl, ...
                     evntMask)
% calculate full duration at specified percentage of amplitude of event
% half maximum value
% t = time
% t_prof = profile
% sE, eE = start and end of event
% peakVal, peakPos, bs = values of peak and baseline
% prcOfAmpl = percentage of amplitude of event where to calculate full
% duration
% eM = mask of event
t_indx = (1:1:numel(t));
t_indx = t_indx(:);

% value of desired percentage of amplitude
prcOfAmpl_val = (peakVal - bs)*(prcOfAmpl/100) + bs;
prcOfAmpl_val_sE = (peakVal - bs)*((min(1,prcOfAmpl-10))/100) + bs;

% parts of profiles before and after maximum of event
t_prof_beforePeak = t_prof(1:peakPos-1);
t_prof_beforePeak = flipud(t_prof_beforePeak(:));

if isempty(t_prof_beforePeak)
    t_prof_beforePeak = t_prof(1);
end

% before peak, also add some points with lower amplitude than desired one
sE = min(sE, ...
    numel(t_prof_beforePeak) - ...
    find((t_prof_beforePeak-prcOfAmpl_val_sE) < 0, 1, 'first'));
% after peak
t_prof_afterPeak = t_prof(peakPos:end);
p_prcOfAmpl_afterPeak = peakPos + ...
    find((t_prof_afterPeak - prcOfAmpl_val) < 0, 1, 'first');
if (isempty(p_prcOfAmpl_afterPeak) || p_prcOfAmpl_afterPeak>numel(t_prof))
    p_prcOfAmpl_afterPeak = numel(t_prof);
end

eE = max(eE, p_prcOfAmpl_afterPeak);

if sE>1
    y_t1 = [ -inf(sE-1,1); t_prof(sE:peakPos-1) ];
    evntMask_1 = [ false(sE-1,1); evntMask(sE:peakPos-1) ];
else
    y_t1 = t_prof(1:peakPos-1);
    evntMask_1 = evntMask(1:peakPos-1);
end
y_t2 = [ -inf(peakPos-1,1); t_prof(peakPos:eE) ];
evntMask_2 = [ false(peakPos-1,1); evntMask(peakPos:eE) ];

% find the first position of desired percentage of amplitude
y_t1_r = flipud(y_t1);
evntMask_1_r = flipud(evntMask_1);
d_y_t1 = gradient(y_t1_r);
if all(d_y_t1==0) || all(isnan(d_y_t1))
    pos_1 = numel(d_y_t1);
else
    indDer1 = find( (d_y_t1<=0 | evntMask_1_r) & y_t1_r<prcOfAmpl_val, ...
        1, 'first')+1;
    if indDer1 > numel(y_t1_r), indDer1 = numel(y_t1_r); end
    if ~isempty(indDer1), y_t1_r(indDer1:end) = -inf; end
    y_t1 = flipud(y_t1_r);
    % get the first position of desired percentage of amplitude
    [~,pos_1] = min(abs(y_t1 - prcOfAmpl_val));
    % pos_last_0_y_t1 = find(y_t1>-inf, 1, 'first')-1;
    % [~,pos_1] = min(abs(y_t1(pos_last_0_y_t1:end) - prcOfAmpl_val));
    % pos_1 = pos_1 + pos_last_0_y_t1 - 1;
end

% find the second position of desired percentage of amplitude
d_y_t2 = gradient(y_t2);
indDer2 = find( (d_y_t2<0 | evntMask_2) & y_t2<prcOfAmpl_val & t_indx(1:eE)>peakPos, ...
    1, 'first')+1;
if indDer2 > numel(y_t2), indDer2 = numel(y_t2); end
if ~isempty(indDer2), y_t2(indDer2:end) = -inf; end

[~, pos_2] = min(abs(y_t2 - prcOfAmpl_val));
% pos_last_numNot0_y_t2 = find(y_t2(peakPos:end)>0, 1, 'last') + peakPos - 1;
% [~, pos_2] = ...
%     min(abs(y_t2(peakPos:pos_last_numNot0_y_t2) - prcOfAmpl_val));
% pos_2 = pos_2 + peakPos - 1;

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