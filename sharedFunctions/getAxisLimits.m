function [limits, best_inc] = getAxisLimits(data, offset, startAt0)
% input: - data 
%        - offset in percentage
%        - start at zero, (true or false)
if ~exist('startAt0', 'var')
    startAt0 = false;
end
data = double(data(:));
if startAt0
    dataRange = abs(max(data)-0);
    min_data = 0;
else
    dataRange = abs(max(data)-min(data));
    min_data = min(data)-dataRange*offset/100;
end
max_data = max(data)+dataRange*offset/100;
% update data range with offset
dataRange = abs(max_data-min_data);
% find nice increment (chosen 5 ticks including min and max)
try
    % get exponent for 5 ticks/increments
    exp_num = floor(log10(dataRange/5));
    % increments
    inc = [1;2;3;5] .* 10.^(exp_num-1:exp_num+1);
    inc = inc(:);
    % find best increment (chosen ~6 ticks including max)
    [~, idx] = min( abs( dataRange ./ inc - 4 ) );
    best_inc = inc(idx);
catch
    best_inc = round(dataRange/3, 1, "significant");
end
% calculate min and max limits
limits = [floor(min_data/best_inc)*best_inc, ...
          ceil(max_data/best_inc)*best_inc];

end





