function params= getParametersOfEventFromRawImage(...
    imgE, imgEs, imgE_m, ...
    bs, t0, r_m, c_m, ...
    t_prof_s, x_prof_s, t, x, imgData, normImgFlag)

% calculate parameters of event from its image

% imgE = event image, not normalized
% imgEs = smoothed event image, 2D spline
% imgE_m = mask of event image
% bs = baseline
% t0 = start of event
% r_m = row where max of event is
% c_m = colum where max of event is
% t_prof_s, x_prof_s = profiles from smoothed image
% t, x = time in ms and spatial coordinate in um of event image

% get parameters of spark from profiles
% half max of image
eventMax = imgEs(r_m,c_m);
if isempty(bs)
    bs = prctile(imgEs(:),5); % take 5th percentile as baseline
    half_max_img = bs + (eventMax-bs)/2;
else
    half_max_img = bs + (eventMax-bs)/2;
end

% half max in profiles
half_max_t = bs + ( t_prof_s(c_m) - bs )/2;
half_max_x = bs + ( x_prof_s(r_m) - bs )/2;

% do treshold of image with half maximum
imgEs_halfAmpl = imgEs;
imgEs_halfAmpl(imgEs < half_max_img) = 0;

% apply mask of detected event, to be sure that we get only parts where
% event was detected
imgEs_halfAmpl(~imgE_m) = 0;

% get properties of detected, connected objects
CC_eventsInSparkRegion = bwconncomp(imgEs_halfAmpl, 8);
% take that one which has maximum determined before inside
% this should be maximum of spark which we want to
% analyze, in case that there are >2 detected
% sparks, which are close to each other and they
% are overlaping into rectangle region of spark
% currently being analyzed
eventsInSparkRegion = regionprops(CC_eventsInSparkRegion,imgEs,...
    'SubarrayIdx','PixelIdxList');
% select the region containing spark maximum
indSpReg = nan;
for r=1:numel(eventsInSparkRegion)
    if ismember(r_m,eventsInSparkRegion(r).SubarrayIdx{1}) && ...
            ismember(c_m,eventsInSparkRegion(r).SubarrayIdx{2})
        indSpReg = r;
    end
end

% if not region found, finish and return nans
if isnan(indSpReg)
    % output
    params.amplitude = nan;
    params.TTP = nan;
    params.FDHM = nan;
    params.FWHM = nan;
    params.sparkMass = nan;
    params.tauD = nan;
    params.bs = nan;
    params.t0 = nan;
    params.t_max = nan;
    params.v_max = nan;

    params.half_max_t = nan;
    params.half_max_t_1 = nan;
    params.half_max_t_2 = nan;
    params.half_max_x = nan;
    params.half_max_x_1 = nan;
    params.half_max_x_2 = nan;
    
    return
end

% region of halfMax event
eventsInSparkRegion = eventsInSparkRegion(indSpReg);

% final mask of halfMax smoothed image
m_imgEs_halfAmpl = false(size(imgEs_halfAmpl));
m_imgEs_halfAmpl(eventsInSparkRegion.PixelIdxList) = true;
% halfMax smoothed image
imgEs_halfAmpl(~m_imgEs_halfAmpl) = 0;

% normalized (F-F0)/(F0-blank)
% take position of maximum as it was determined
% before, take amplitude from time profile
if ~normImgFlag
    eventAmpl = ( t_prof_s(c_m)-bs ) / ( bs-imgData.blank );
else
    eventAmpl = t_prof_s(c_m) - bs;
end

% find t0, if does not exist
if isempty(t0)
    t0_ind_est = find( imgEs(r_m,1:c_m)<bs,1,'last');
    if isempty(t0_ind_est) || t0_ind_est<1
        t0_ind_est = eventsInSparkRegion.SubarrayIdx{2}(1);
    end
    t0 = t(t0_ind_est);    
end

if t(c_m) < t0
    t0 = t(c_m);
end

TTP = t(c_m)-t0;

% use profiles from tresholded image
pT_halfMax = imgEs_halfAmpl(r_m,:)>0;
pX_halfMax = imgEs_halfAmpl(:,c_m)>0;

t_prof_halfMax = t_prof_s(pT_halfMax);
x_prof_halfMax = x_prof_s(pX_halfMax);

% positions of half max in t axis
[~,half_max_t_1] = min( abs(t_prof_s(eventsInSparkRegion.SubarrayIdx{2}(1):c_m) - half_max_t) );
half_max_t_1 = half_max_t_1 + eventsInSparkRegion.SubarrayIdx{2}(1)-1;
if isempty(half_max_t_1), half_max_t_1=1; end

[~,half_max_t_2] = min( abs(t_prof_s(c_m:eventsInSparkRegion.SubarrayIdx{2}(end)) - half_max_t) );
half_max_t_2 = half_max_t_2 + c_m -1;
if isempty(half_max_t_2), half_max_t_2=numel(t_prof_halfMax); end

half_max_t_1 = t(half_max_t_1);
half_max_t_2 = t(half_max_t_2);

FDHM = half_max_t_2 - half_max_t_1;

% positions of half max in x axis
[~,half_max_x_1] = min( abs(x_prof_s(eventsInSparkRegion.SubarrayIdx{1}(1):r_m) - half_max_x) );
half_max_x_1 = half_max_x_1 + eventsInSparkRegion.SubarrayIdx{1}(1)-1;
if isempty(half_max_x_1), half_max_x_1=1; end

[~,half_max_x_2] = min( abs(x_prof_s(r_m:eventsInSparkRegion.SubarrayIdx{1}(end)) - half_max_x) );
half_max_x_2 = half_max_x_2 + r_m -1;
if isempty(half_max_x_2), half_max_x_2=numel(x_prof_halfMax); end

half_max_x_1 = x(half_max_x_1);
half_max_x_2 = x(half_max_x_2);

FWHM = half_max_x_2 - half_max_x_1;

sparkMass = eventAmpl*1.206*FWHM^3;

% estimate tauD
profD = t_prof_s(c_m+1:end);
try
    bs_end = mean(profD(end-2:end)); % last 3 points
catch
    bs_end = profD(end); % last point
end
[~,p_tauD] = min( abs( profD - (0.33*( imgEs(r_m,c_m) - bs_end)+bs_end) ) );
tauD = p_tauD * imgData.pxSzT;

% output
params.amplitude = eventAmpl;
params.TTP = TTP;
params.FDHM = FDHM;
params.FWHM = FWHM;
params.sparkMass = sparkMass;
params.tauD = tauD;
params.bs = bs;
params.t0 = t0;
params.t_max = t(c_m);
params.v_max = t_prof_s(c_m);

params.half_max_t = half_max_t;
params.half_max_t_1 = half_max_t_1;
params.half_max_t_2 = half_max_t_2;
params.half_max_x = half_max_x;
params.half_max_x_1 = half_max_x_1;
params.half_max_x_2 = half_max_x_2;

end

