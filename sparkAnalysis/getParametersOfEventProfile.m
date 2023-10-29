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
        % half maximum value
        half_max_t = (val_t - bs_t)*0.5 + bs_t;
        max_t_25perc = (val_t - bs_t)*0.25 + bs_t;
        max_t_75perc = (val_t - bs_t)*0.75 + bs_t;

        % parts of profiles before and after maximum of event
        % extend, to be sure points with 25% and 50% of max values are there
        t_prof_beforePeak = t_prof(1:pos_t-1);
        t_prof_beforePeak = flipud(t_prof_beforePeak(:));
        % before peak
        sE = min(sE, ...
            numel(t_prof_beforePeak) - ...
            find((t_prof_beforePeak-max_t_25perc) < 0, 1, 'first'));
        % after peak
        t_prof_afterPeak = t_prof(pos_t:end);
        p_50_afterPeak = pos_t + ...
            find((t_prof_afterPeak-half_max_t) < 0, 1, 'first');
        if (isempty(p_50_afterPeak) || p_50_afterPeak>numel(t_prof))
            p_50_afterPeak = numel(t_prof);
        end
        eE = max(eE, p_50_afterPeak);

        if sE>1
            y_t1 = [ zeros(sE-1,1); t_prof(sE:pos_t-1) ];
        else
            y_t1 = t_prof(1:pos_t-1);
        end
        y_t2 = [ zeros(pos_t-1,1); t_prof(pos_t:eE) ];

        % find position of first half max position
        y_t1_r = flipud(y_t1);
        d_y_t1 = gradient(y_t1_r);

        if all(d_y_t1==0)
            pos_25 = numel(d_y_t1);
            pos_75 = numel(d_y_t1);
            pos_50_1 = numel(d_y_t1);
        else
            indDer1_25 = find( d_y_t1<=0 & y_t1_r<max_t_25perc, 1, 'first');
            indDer1_50 = find( d_y_t1<=0 & y_t1_r<half_max_t, 1, 'first');

            y_t1_r_25 = y_t1_r;
            y_t1_r_50 = y_t1_r;
            if ~isempty(indDer1_25), y_t1_r_25(indDer1_25:end) = 0; end
            if ~isempty(indDer1_50), y_t1_r_50(indDer1_50:end) = 0; end

            y_t1_25 = flipud(y_t1_r_25);
            y_t1_50 = flipud(y_t1_r_50);

            [~,pos_50_1] = min(abs(y_t1_50 - half_max_t));
            % find also positions of 25 and 75% of amplitude
            [~,pos_25] = min(abs(y_t1_25 - max_t_25perc));
            [~,pos_75] = min(abs(y_t1_50 - max_t_75perc));
        end

        % find position of second half max position
        d_y_t2 = gradient(y_t2);

        indDer2 = find( d_y_t2<0 & y_t2<(max_t_25perc) & t_indx(1:eE)>pos_t, ...
            1, 'first');
        if ~isempty(indDer2)
            y_t2(indDer2:end) = 0;
        end

        [~,pos_50_2] = min(abs(y_t2 - half_max_t));

        % positions in time
        half_max_t_1 = t(pos_50_1);
        half_max_t_2 = t(pos_50_2);
        % get full duration in half maximum
        FDHM = half_max_t_2 - half_max_t_1;
        if FDHM<0
            FDHM=0;
        elseif isempty(FDHM)
            FDHM=0;
        end

    catch
        FDHM = w_t(1);
        half_max_t = nan;
        half_max_t_1 = nan;
        half_max_t_2 = nan;
        pos_25 = nan;
        pos_75 = nan;
    end

    % find t0 (start of spark) as a cross section of 0 and line fit of part
    % of event with amplitudes between 25 and 75 % of max amplitude

    if pos_25 == pos_75
        pos_25 = pos_25 - 1;
    end

    if isempty(t0)
        try
            f_line25_75 = fit(t(pos_25:pos_75),t_prof(pos_25:pos_75),'poly1');
            t_fit = t(find(t_prof>bs_t,1,'first'):pos_t);
            line25_75 = feval(f_line25_75,t_fit);
            t0 = t_fit(find(line25_75<bs_t,1,'last'));
        catch
            % find t0 (start of spark) as 5% of   increase of  set manually
            % t0 = t( find(t_prof > ((val_t-bs_t)/50 + bs_t), 1, 'first'));
            t0 = t( numel(t_prof) - ...
                find(flipud(t_prof) < ((val_t-bs_t)/50 + bs_t), 1, 'first') + 1);
        end

        %     if isempty(t0)
        %         % find t0, start of spark, 0.05 set manually, not really nice
        %         t0 = t( find(t_prof > ((val_t-bs_t)/50 + bs_t),1,'first'));
        %     end

    end

    % get time to peak
    TTP = t(pos_t) - t0;
    if TTP<0
        TTP=0;
    elseif isempty(TTP)
        TTP=0;
    end
else
    FDHM = nan;
    half_max_t = nan;
    half_max_t_1 = nan;
    half_max_t_2 = nan;
    pos_25 = nan;
    pos_75 = nan;
    t0 = nan;
    TTP = nan;
end


%% get spatial parameters
if ~isempty(x_prof)
    %[val_x,pos_x] = max(x_prof);
    [pks_x,locs_x,w_x,~] = findpeaks(x_prof, x, ...
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
        half_max_x = (val_x(1)-bs_x)/2 + bs_x;
        max_x_25perc = (val_x(1)-bs_x)*0.25 + bs_x;

        y_x1 = x_prof(1:pos_x-1);
        y_x2 = cat(1 ,zeros(pos_x-1,1), x_prof(pos_x:(length(x_prof))));

        % find position of first half max
        y_x1_r = flipud(y_x1);
        d_y_x1 = gradient(y_x1_r);

        indDer1 = find( d_y_x1>=0 & y_x1_r<half_max_x,1,'first');
        if ~isempty(indDer1)
            y_x1_r(indDer1:end) = 0;
        end

        y_x1 = flipud(y_x1_r);

        [~,pos1_x] = min(abs(y_x1 - half_max_x));

        % find position of second half max
        d_y_x2 = gradient(y_x2);

        indDer2 = find( d_y_x2>=0 & y_x2<max_x_25perc & x_indx>pos_x,1,'first');
        if ~isempty(indDer2)
            y_x2(indDer2:end) = 0;
        end

        [~,pos2_x] = min(abs(y_x2 - half_max_x));

        half_max_x_1 = x(pos1_x);
        half_max_x_2 = x(pos2_x);

        FWHM = half_max_x_2 - half_max_x_1;

        if isempty(FWHM) || isnan(FWHM)
            FWHM = w_x(1);
        end

        if FWHM > max(x), FWHM = max(x); end

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
eventParams.FDHM = FDHM;
eventParams.FWHM = FWHM;
eventParams.sparkMass = sparkMass;
eventParams.bs_t = bs_t;
eventParams.bs_x = bs_x;
eventParams.t0 = t0;
eventParams.t_max = t(pos_t);
eventParams.v_max = val_t;
eventParams.half_max_t = half_max_t;
eventParams.half_max_t_1 = half_max_t_1;
eventParams.half_max_t_2 = half_max_t_2;
eventParams.half_max_x = half_max_x;
eventParams.half_max_x_1 = half_max_x_1;
eventParams.half_max_x_2 = half_max_x_2;

end