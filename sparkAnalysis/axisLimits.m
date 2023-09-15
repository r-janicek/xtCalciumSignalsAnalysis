function limits = axisLimits(x)
% input is a vector
% find reasonable minimum and maximum for axis limits

% increments
inc = [0.01,0.03,0.05, 0.1,0.2,0.3,0.5, 1,2,3,5, 10,20,50, ...
    100,200,300,500, 1000,2000,3000,5000, 10000,30000,50000,100000];

% find best increment (chosen 7 ticks including min and max)
[~, idx] = min( abs( (max(x)-min(x)) ./ inc - 5 ) );

% find min and max limits
max_lim = ceil(max(x)/inc(idx))*inc(idx);
min_lim = floor(min(x)/inc(idx))*inc(idx);

limits = [min_lim, max_lim];
 
end