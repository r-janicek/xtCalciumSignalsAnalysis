function analysisOpt(hO,~,mainFig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

hObjs = getappdata(mainFig,'hObjs');

switch hO.String

    case 'show individual events'
        
        if hO.Value
            hObjs.check_sparkParams.Value = 1;                   
        end
        
    case 'get sparks params'
        
        if ~hO.Value
            hObjs.check_showEventsFigs.Value = 0;                   
        end
        
end


end

