function limits = getAxisLimits(data, offset)
    % data
    % offset in percentage
    data = data(:);
    dataRange = abs(max(data)-min(data));
    limits = [min(data)-dataRange*offset/100, ...
        max(data)+dataRange*offset/100];
end