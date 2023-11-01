function limits = getAxisLimits(data, offset, startAt0)
% input: - data 
%        - offset in percentage
%        - start at zero, (true or false)
if ~exist('startAt0', 'var')
    startAt0 = false;
end
data = data(:);
if startAt0
    dataRange = abs(max(data)-0);
    min_data = 0;
else
    dataRange = abs(max(data)-min(data));
    min_data = min(data)-dataRange*offset/100;
end
max_data = max(data)+dataRange*offset/100;
% find nice increment (chosen 5 ticks including min and max)
inc = round(dataRange/3, 1, "significant");
% calculate min and max limits
limits = [floor(min_data/inc)*inc, ...
          ceil(max_data/inc)*inc];
end






