function eventParams = getParametersOfEventProfile( ... 
    t, t_prof, x, x_prof, ...
    bs_t, bs_x, t0, blank, peakData, tProf_m, normImgFlag)

% bs = baseline
% t0 = start of spark
% blank 

if isempty(tProf_m)
    tProf_m = true(size(t));
end

t_indx = (1:1:length(t));
t_indx = t_indx(:);

x_indx = (1:1:length(x));
x_indx = x_indx(:);

% start and end of event from event t profile mask
sE = find(tProf_m,1,'first');
if isempty(sE), sE = 1; end
eE = find(tProf_m,1,'last');
if isempty(eE), eE = numel(tProf_m); end

% (F-F0)/(F0-blank)
if normImgFlag
    if isempty(peakData)
        Ampl = max(t_prof)-bs_t;
    else
        Ampl = peakData.val-bs_t;
    end
else
    if isempty(peakData)
        Ampl = (max(t_prof)-bs_t)/(bs_t-blank);
    else
        Ampl = (peakData.val-bs_t)/(bs_t-blank);
    end
end

%% time parameters
if ~isempty(t_prof)
    % search for peak
    if isempty(peakData)
        %[val_t,pos_t] = max(t_prof);
        [pks_t,locs_t,w_t,~] = findpeaks(t_prof(tProf_m),t(tProf_m),...
            'SortStr','descend');
    else
        % look 15 ms around, in points
        eps = ceil(15/mean(diff(t)));
        startPoint = peakData.pos-eps;
        if startPoint<1, startPoint=1; end
        endPoint = peakData.pos+eps;
        if endPoint>numel(t_prof), endPoint=numel(t_prof); end

        [pks_t,locs_t,w_t,~] = findpeaks( ...
            t_prof(startPoint:endPoint),...
            t(startPoint:endPoint),'SortStr','descend');
    end
    try
        val_t = pks_t(1);
        [~, pos_t] = min(abs(t-locs_t(1)));
    catch
        pks_t = t_prof(1);
        locs_t = 1;
        w_t = 0;
        val_t = t_prof(1);
        pos_t = 1;
    end

    try
        % get full duration in half maximum
        [FDHM, half_max_t, half_max_t_1, half_max_t_2] = ...
            fullDurationCalc(t, t_prof, sE, eE, val_t, pos_t, bs_t, 50, tProf_m);
        % get positions at defined percentiles (part of event rise fitted by line)
        percRise = [10, 90];
        [~, ~, ~, ~, pos_percRise_1, ~] = fullDurationCalc( ...
            t, t_prof, sE, eE, val_t, pos_t, bs_t, percRise(1), tProf_m);
        [~, ~, ~, ~, pos_percRise_2, ~] = fullDurationCalc( ...
            t, t_prof, sE, eE, val_t, pos_t, bs_t, percRise(2), tProf_m);
        % [~, ~, ~, ~, pos_10, ~] = fullDurationCalc( ...
        %     t, t_prof, sE, eE, val_t, pos_t, bs_t, 10);
        % [~, ~, ~, ~, pos_90, ~] = fullDurationCalc( ...
        %     t, t_prof, sE, eE, val_t, pos_t, bs_t, 90);
    catch
        FDHM = w_t(1);
        half_max_t = nan;
        half_max_t_1 = nan;
        half_max_t_2 = nan;
        pos_percRise_1 = nan;
        pos_percRise_2 = nan;
    end

    % find t0 (start of spark) as a cross section of 0 and line fit of part
    % of event with amplitudes between 25 and 75 % of max amplitude
    if pos_percRise_1 == pos_percRise_2
        pos_percRise_1 = pos_percRise_1 - 1;
    end
    if pos_percRise_1 == pos_percRise_2
        pos_percRise_1 = pos_percRise_1 - 1;
    end

    try
        f_line_percRise = fit( ...
            t(pos_percRise_1:pos_percRise_2), ...
            t_prof(pos_percRise_1:pos_percRise_2), 'poly1');
        t0_fittedLine = (bs_t-f_line_percRise.p2)/f_line_percRise.p1;
        assert(t0_fittedLine>0)
    catch
        % find t0 (start of spark) as 1% amplitude increase over baseline
        t0_fittedLine = t( numel(t_prof) - ...
            find(flipud(t_prof) < ((val_t-bs_t)/100 + bs_t), 1, 'first') + 1);
    end
    if isempty(t0)
        t0 = t0_fittedLine; 
    end
    % get time to peak
    TTP = t(pos_t) - t0;
    TTP_fittedLine = t(pos_t) - t0_fittedLine;
    if (TTP<0) || isempty(TTP), TTP=0; end
    if (TTP_fittedLine<0) || isempty(TTP_fittedLine), TTP_fittedLine=0; end
    
else
    FDHM = nan;
    half_max_t = nan;
    half_max_t_1 = nan;
    half_max_t_2 = nan;
    pos_percRise_1 = nan;
    pos_percRise_2 = nan;
    t0 = nan;
    t0_fittedLine = nan;
    TTP = nan;
    TTP_fittedLine = nan;
end


%% get spatial parameters
if ~isempty(x_prof)
    %[val_x,pos_x] = max(x_prof);
    [pks_x, locs_x, w_x, ~] = findpeaks(x_prof, x, ...
        'SortStr','descend', 'Annotate','extents');
    try
        if pks_x == max(x_prof)
            val_x = pks_x(1);
            [~, pos_x] = min(abs(x-locs_x(1)));
            %pos_x = round(locs_x(1)/mean(diff(x)));
        else
            [val_x, pos_x] = max(x_prof);
            w_x = round(numel(x_prof)/2)*mean(diff(x));
        end
    catch
        pks_x = x_prof(1);
        locs_x = 1;
        w_x = 0;
        val_x = pks_x(1);
        pos_x = 1;
    end

    try
        % get full width in half maximum
        [FWHM, half_max_x, half_max_x_1, half_max_x_2] = ...
            fullDurationCalc(x, x_prof, 1, numel(x_prof), ...
            val_x, pos_x, bs_x, 50, true(size(x)));
    catch
        FWHM = w_x(1);
        half_max_x = nan;
        half_max_x_1 = nan;
        half_max_x_2 = nan;
    end
else
    FWHM = nan;
    half_max_x = nan;
    half_max_x_1 = nan;
    half_max_x_2 = nan;
end

%calculating spark mass
sparkMass = Ampl*1.206*FWHM^3;

eventParams.amplitude = Ampl;
eventParams.TTP = TTP;
eventParams.(matlab.lang.makeValidName( ...
    sprintf('TTP_lineFit_%d_%d',percRise(1),percRise(2)))) = TTP_fittedLine;
%eventParams.TTP_fittedLine = TTP_fittedLine;
eventParams.FDHM = FDHM;
eventParams.FWHM = FWHM;
eventParams.sparkMass = sparkMass;
eventParams.bs_t = bs_t;
eventParams.bs_x = bs_x;
eventParams.t0 = t0;
eventParams.(matlab.lang.makeValidName( ...
    sprintf('t0_lineFit_%d_%d',percRise(1),percRise(2)))) = t0_fittedLine;
%eventParams.t0_fittedLine = t0_fittedLine;
eventParams.t_max = t(pos_t);
eventParams.v_max = val_t;
eventParams.half_max_t = half_max_t;
eventParams.half_max_t_1 = half_max_t_1;
eventParams.half_max_t_2 = half_max_t_2;
eventParams.half_max_x = half_max_x;
eventParams.half_max_x_1 = half_max_x_1;
eventParams.half_max_x_2 = half_max_x_2;

end