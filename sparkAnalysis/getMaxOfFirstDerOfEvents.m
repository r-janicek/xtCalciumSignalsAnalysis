function [firstDerMaxValPos,firstDerMaxVal] = getMaxOfFirstDerOfEvents(prof,startOfEvent,posOfPeaks)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% 1.derivative
prof_Der = gradient(prof);

firstDerMaxValPos = zeros(size(startOfEvent));
firstDerMaxVal = zeros(size(startOfEvent));
for i=1:numel(posOfPeaks)
    
    try
        % part of derivative
        prof_Der_p = prof_Der(startOfEvent(i):posOfPeaks(i));
        
        [v,p] = max(prof_Der_p);
        if p < round(numel(prof_Der_p)/2)
            [v,p] = max([zeros(round(numel(prof_Der_p)/2)-1,1) ;prof_Der_p(round(numel(prof_Der_p)/2):end)]);
        end
    catch
        v = nan;
        p = nan;
    end
    
    firstDerMaxValPos(i) = startOfEvent(i) + p - 1;
    firstDerMaxVal(i) = v;
    
    clearvars prof_Der_p prof_p v p
    
end

end

