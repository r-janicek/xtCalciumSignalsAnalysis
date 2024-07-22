function analysisOpt(hO,~,mainFig)
hObjs = getappdata(mainFig,'hObjs');
switch hO.String
    case 'show individual events'
        if hO.Value
            hObjs.check_sparkParams.Value = 1;
            if isfield(getappdata(mainFig), 'sparkDetection')
                sparkDetection = getappdata(mainFig, 'sparkDetection');
                analyzedEvntsBrowserTbl = ...
                    sparkDetection.analyzedEvntsBrowserTbl;
                if ~isempty(analyzedEvntsBrowserTbl)
                    % close previous events viewer window, if any
                    close(findall(0, 'Type','Figure', 'Tag','CaEventsBrowser'))
                    eventsBrowserWindow(mainFig, analyzedEvntsBrowserTbl)
                end
            end
        end
        
    case 'get sparks params'
        if ~hO.Value
            hObjs.check_showEventsFigs.Value = 0;                   
        end 
end
end

