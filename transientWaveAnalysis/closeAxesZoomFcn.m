function closeAxesZoomFcn(hO,E,analysisFig)

% delete big axes figure
delete(hO)

% do fitting and analyzing of profile, in small axes in analysis window
findAndAnalyzePeaks([],[],analysisFig)

end

