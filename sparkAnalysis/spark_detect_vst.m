function [bw, im_detect] = spark_detect_vst(im, filt_space_fwhm, smooth_iter, base_iter, do_vst, tau, expansion_factor)
% 
% Calcium spark detection in linescan images using variance stabilization and local baseline subtraction.
% The code has been written using MATLAB R2010a, 
% but can also be run without modification within Octave (open source, available to download at http://www.gnu.org/software/octave/).
% It has been tested using Octave version 3.2.4.
%
% 
% 
% Input parameters:
% IM - A 2D matrix representing the pixel data. Any constant offset
% should have been subtracted in advance.
% IMPORTANT: The time dimension is assumed to be horizontal, the spatial dimension is vertical.
% FILT_SPACE_FWHM - The full-width-at-half-maximum (FWHM) of the Gaussian 
%                   spatial filter. This should be <= the FWHM of a typical spark, defined in terms of pixels. (Default = 8)
% SMOOTH_ITER - The number of smoothing iterations to apply along the
%               temporal dimension. Higher values suppress noise, but excessively
%               high values cause sparks to be suppressed too. (Default = 3) 
% BASE_ITER - The number of smoothing iterations to apply for estimating
%             the baseline. The must be greater than SMOOTH_ITER. 
%             If sparks occur from a constant baseline, high values can be used 
%             (effectively resulting in the mean fluorescence of the row being used to estimate the baseline). 
%             Because the filter size increases dramatically, this shouldn't be higher than about 10. 
%             For quickly varying baselines, smaller values are required to avoid undesirable detections. (Default = 5)
% DO_VST - TRUE if the square root transform should be applied for variance stabilization, FALSE otherwise. 
%          If the baseline varies a lot, TRUE should always be chosen. 
%          For constant baselines, FALSE can give more detections, but one is less confident that they are correct. (Default = TRUE)
% TAU - Threshold value, defined in terms of noise standard deviations. (Default = 5)
% EXPANSION_FACTOR - A value between 0 and 1 that determines the size of detected spark regions. 
%                    Lower values give larger regions. 
%                    A value of 1 means that sparks are detected by simply thresholding the processed image at TAU. 
%                    Such regions might be too small to contain the peaks of events. (Default = 0.5)
% 
%
% Output parameters:
% BW - A binary image, the same size as IM, containing detected spark regions.
% IM_DETECT - The processed, baseline-subtracted image that is thresholded to give BW. 
%             Pixel values in this are normalized to the estimated noise, 
%             and so correspond to the threshold that would be required to detect that particular pixel as significant.
% 
% 
% Dependencies for MATLAB:
% The function IMFILTER from the Image Processing Toolbox is required. 
% If this toolbox is unavailable, the dependency can be removed by replacing calls to IMFILTER with CONV2. 
% However, note that in at least some MATLAB versions CONV2 is much less efficient for large filters with many zero elements,
% so the advantages of the iterative algorithm are lost. 
% Also, CONV2 uses zero-padding - so reduces the usefulness of results at image boundaries.
% 
% Dependencies for Octave:
% Much of the Image Processing Toolbox, including IMFILTER, is replicated in the IMAGE package at http://octave.sourceforge.net/image/
% Running the code:
% The contents of this file should be copied to a text file named 'spark_detect_vst.m',
% which is added to the search path for MATLAB or Octave.
% It can then be run, for example, by typing the following commands: im = imread('test_image.tif'); 
% Read an original linescan image imagesc(im); % Display the image
% [bw, im_detect] = spark_detect_vst(im, 10, 3, 6); Apply detection 
% figure; imagesc(bw); % Show detected regions
% 
% AUTHOR: Peter Bankhead, while at Queen's University, Belfast.
% Email contact: peter.bankhead{at}bioquant{dot}uni-heidelberg{dot}de
% COPYRIGHT: The source code and algorithm are to be used for non-commercial purposes only. For commercial use, contact the author.



%% CHECK INPUT ARGUMENTS --------------------------------------------------------------------------

% Ensure first argument is an image
if nargin < 1 || ~isnumeric(im) || isvector(im) || ndims(im) ~= 2
    error('A valid 2D image was not passed to SPARK_DETECT_VST!');
else
    % Make sure the image is floating point
    if ~isfloat(im)
        im = single(im);
    end
end

% Check Gaussian filter FWHM
if nargin < 2 || isempty(filt_space_fwhm)
    filt_space_fwhm = 8;
elseif ~isscalar(filt_space_fwhm) || filt_space_fwhm < 0
    error('FILT_SPACE_FWHM must a positive scalar');
end

% Check smoothing iterations
if nargin < 3 || isempty(smooth_iter)
    smooth_iter = 3;
elseif ~isscalar(smooth_iter) || smooth_iter < 0
    error('SMOOTH_ITER must be a positive integer');
end

% Check baseline smoothing iterations
if nargin < 4 || isempty(base_iter)
    base_iter = 5;
elseif ~isscalar(base_iter) || base_iter <= smooth_iter
    error('BASE_ITER must be an integer greater than SMOOTH_ITER');
end

% Do variance stabilization by default
if nargin < 5 || isempty(do_vst)
    do_vst = true;
elseif ~isscalar(do_vst) && ~islogical(do_vst)
    error('DO_VST should be TRUE or FALSE');
end

% Check threshold
if nargin < 6 || isempty(tau)
    tau = 5;
elseif ~isscalar(tau) && ~isnumeric(tau) || tau <= 0
    error('TAU must be a positive numeric scalar');
end

% Check the region expansion factor
if nargin < 7
    expansion_factor = 0.5;
elseif ~isscalar(expansion_factor) && ~isnumeric(expansion_factor) && ...
        expansion_factor < 0 || expansion_factor > 1 
    error('EXPANSION_FACTOR must be a numeric scalar between 0 and 1');
end



%% APPLY SPATIAL FILTERING AND VARIANCE STABILIZATION----------------------

% Generate a 1D Gaussian filter
filt_space = gaussian_filt_1d(filt_space_fwhm);

% Apply spatial filtering
im = imfilter(im, filt_space(:), 0);

% Apply the square root transform for variance stabilization (if desired) 
if do_vst
    im = realsqrt(im);
end



%% ESTIMATE NOISE----------------------------------------------------------

% Estimate the noise in the filtered and variance-stabilized image.
% This shouldn't include pixels that included filtering over boundaries. 
b_pad = ceil(numel(filt_space)/2);
noise_raw_std = mad_noise_estimate(im(b_pad:end-b_pad+1,:));



%% APPLY TEMPORAL SMOOTHING AND BASELINE ESTIMATION------------------------

% Use a filter derived from the cubic B-spline by default 
b3 = [1, 4, 6, 4, 1] / 16;

% Initialize the direct smoothing filter 
filt_direct = 1;

% Iterate through the smoothing, upsampling the filter by inserting
% 2^(ii-1)-1 zeros between filter coefficients at each iteration - as in a
% stationary wavelet transform.
for ii = 1:base_iter
    % Apply the filtering
    filt_temp = insert_filter_zeros(b3, 2^(ii-1)-1);
    im = imfilter(im, filt_temp, 'symmetric');
    
    % Update the direct filter - this could be applied to the image
    % directly to get the same result instead of using the (more efficient)
    % iterative process
    filt_direct = conv(filt_direct, filt_temp);
    
    % If the initial smoothing is complete, keep this image and direct filter
    if ii == smooth_iter
        im_smooth = im;
        filt_time_smooth = filt_direct;
    end
    
end

% The final smoothed result is the baseline image 
im_base = im;
filt_time_base = filt_direct;



%% SUBTRACT THE BASELINE AND UPDATE THE NOISE ESTIMATE---------------------

% Compute the baseline-subtracted detection image 
im_detect = im_smooth - im_base;
% Compute the 'direct' filter for the subtraction image after padding the 
% filters to be the same size
d = numel(filt_time_base) - numel(filt_time_smooth);
filt_time_smooth = [zeros(1, d/2), filt_time_smooth, zeros(1, d/2)]; 
filt_time = filt_time_smooth - filt_time_base;

% Update noise estimate according to the direct filter 
noise_std = noise_raw_std * sqrt(sum(filt_time(:).^2));



%% NORMALIZE DETECTION IMAGE AND THRESHOLD---------------------------------

% Normalize the detection image to the estimated noise 
im_detect = im_detect / noise_std;
% Apply the threshold
bw = region_growing_threshold(im_detect, tau, expansion_factor);



%% ------------------------ SUBFUNCTIONS ----------------------------------
%--------------------------------------------------------------------------

%% ADAPTIVE THRESHOLD------------------------------------------------------ 
function bw = region_growing_threshold(im_detect, tau, expansion_factor)
% Locate local maxima with values MAX_VAL > T_HIGH, and expand these until 
% surrounded by pixels less than MIN(EXPANSION_FACTOR * MAX_VAL, T_HIGH).
% In other words, region sizes depend upon a local threshold defined as a
% factor of the region's maximum value, but which can't be any higher than
% T_HIGH.
% Multiple local maxima can end up in the same region. Therefore, the
% algorithm is applied to each maximum in descending order until all
% maxima exceeding the threshold have been assigned to a region.
% Pad the image with values that definitely won't exceed the threshold -
% this makes it a bit easier when getting the neighbours of above-threshold 
% peaks
im_detect = im_detect([1 1:end end], [1 1:end end]);
im_detect([1, end], :) = -inf;
im_detect(:, [1, end]) = -inf;

% If the Image Processing Toolbox function PADARRAY is available, this would be
% im_detect = padarray(im_detect, [1, 1], -Inf);

% If a pixel has a value >= 2*T_HIGH, then the lower threshold will always
% be T_HIGH, rather than half the amplitude - so IMRECONSTRUCT could be 
% used for a simple hysteresis threshold to process these regions first.
% However, IMRECONSTRUCT isn't currently implemented in Octave, so for 
% compatibility we don't use it.
% Still, the appropriate commands would be this:
% bw = imreconstruct(im_detect >= tau/expansion_factor, im_detect >= tau); 
% Remove any detected regions from IM_DETECT
% im_detect(bw) = -Inf;
bw = false(size(im_detect));

% Compute the offsets added to each linear index to get the pixel's (8) neighbors
m = size(im_detect, 1);
inds_offsets = [-m-1, -m, -m+1, -1, 1, m-1, m, m+1]; 
% Find and sort all above-threshold pixels
inds_all = find(im_detect(:) >= tau);
[dum, s_inds] = sort(im_detect(inds_all), 'descend'); 
inds_all = inds_all(s_inds);

% While there is still anything left in INDS_ALL, expand the next peak
while ~isempty(inds_all)
    % Compute the threshold for this particular region (as a proportion of
    % the peak value). The MIN command wouldn't be necessary if
    % IMRECONSTRUCT was used earlier, since all peaks being processed here
    % would be small enough.
    inds = inds_all(1);
    thresh = min(im_detect(inds) * expansion_factor, tau);
    
    % While we've still got indices for above-threshold pixels, process
    while ~isempty(inds)
        % Update the binary image to include above-threshold pixels
        bw(inds) = true;
        % Remove from detection image
        im_detect(inds) = -Inf;
        % Check the 8-neighbors for being above-threshold 
        inds_temp = bsxfun(@plus, inds(:), inds_offsets); 
        inds_temp(im_detect(inds_temp) < thresh) = [];
        % Update the array of indices for above-threshold pixels only 
        inds = unique(inds_temp(:));
    end
    
    % Remove all the maxima from INDS_ALL that have now been assigned
    inds_all(bw(inds_all)) = [];

end

% Remove padding from the binary image 
bw = bw(2:end-1, 2:end-1);



%% CREATE A 1D GAUSSIAN FILTER---------------------------------------------
function g = gaussian_filt_1d(fwhm)
% Create a 1D Gaussian filter with full-width-at-half-maximum FWHM. 
% Coefficients are normalized so that their sum is 1.

% Convert FWHM to Gaussian sigma value 
sigma = fwhm ./ sqrt(8 * log(2));

% Choose a suitable length - should be an odd number 
len = floor(sigma * 3) * 2 + 1;

% Find where to evaluate the Gaussian 
xx = 0:len-1;
xx = xx - xx(end)/2;

% Evaluate the filter coefficients 
g = exp(-xx.^2 / (2 * sigma.^2));

% Normalize the filter 
g = g / sum(g(:));



%% COMPUTE HAAR MEDIAN ABSOLUTE DEVIATION NOISE ESTIMATE-------------------
function sigma = mad_noise_estimate(im)
% Compute a noise estimate based upon the Median Absolute Deviation of Haar 
% wavelet coefficients of IM, computed along the temporal dimension.
% (Note: A more optimized median computation would help here for very large 
% images, otherwise a subset of the pixels may be used.)

% Compute differences along the time dimension 
diffs = diff(im, [], 2);
% Compute absolute values of required coefficients 
abs_diffs = abs(diffs(:));
% Compute noise estimate
sigma = median(abs_diffs) / sqrt(2) / 0.6745;



%% INSERT ZEROS BETWEEN FILTER COEFFICIENTS--------------------------------
function h2 = insert_filter_zeros(h, n_zeros)
% Insert N_ZEROS zeros between the coefficients of a filter H. 
% H2 is the new, expanded filter.

% Determine the length of the new filter 
len = n_zeros * (numel(h)-1) + numel(h); 
spacing = n_zeros+1;

% Put in the filter coefficients 
h2 = zeros(1, len); 
h2(1:spacing:len) = h;


