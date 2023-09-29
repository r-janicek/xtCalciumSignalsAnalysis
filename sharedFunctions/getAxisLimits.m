function limits = getAxisLimits(data, offset)
% input: data and offset in percentage
data = data(:);
dataRange = abs(max(data)-min(data));
min_data = min(data)-dataRange*offset/100;
max_data = max(data)+dataRange*offset/100;
% find nice increment (chosen 7 ticks including min and max)
inc = round(dataRange/5, 1, "significant");
% calculate min and max limits
limits = [floor(min_data/inc)*inc, ...
          ceil(max_data/inc)*inc];

end