function expDFitSelectionButtonDownFcn(hO,E,analysisFig)

try   
    hO.Parent.Parent.WindowButtonMotionFcn = {@wbmcb,hO.Parent,hO};
    hO.Parent.Parent.WindowButtonUpFcn = {@wbucb,hO.Parent,hO,analysisFig};
    
catch
    return
end
 

% function executed when move with mouse pointer
function wbmcb(hf,~,ha,hO)
        
   cp = ha.CurrentPoint;
   
   % set new position of object
   hO.XData = cp(1,1);
   hO.YData = cp(1,2);
        
   drawnow      
end

% stop function, when mouse button is released
function wbucb(fig,~,ha,hO,analysisFig)
    
fig.WindowButtonMotionFcn = '';
fig.WindowButtonUpFcn = '';

% find position of marker in x and y
hl_splFit = findall(ha,'Type','line','Tag','splineFit');
[~,pos] = min( abs(hl_splFit.XData - hO.XData) );

% in case of peak, search for max in +- 100 ms
if strcmp(hO.Tag,'peakFitSelPoints')

    hl_splFit = findall(ha,'Type','line','Tag','splineFit');

    % look +- 100 ms and find peak if any
    r = ceil(100/mean(diff(hl_splFit.XData)));
    try
        [p_v,p_p] = findpeaks(hl_splFit.YData(pos-r:pos+r),'SortStr','descend');
        p_v = p_v(1);
        p_p = p_p(1);
    catch
        [p_v,p_p] = max(hl_splFit.YData(pos-r:pos+r));
    end
    
    % set new position of object
    hO.XData = hl_splFit.XData(pos - r + p_p - 2);
    hO.YData = p_v;
    
else
    % update position of marker
    hO.XData = hl_splFit.XData(pos);
    hO.YData = hl_splFit.YData(pos);
end
       
% save data
hBsSelCirc = findall(ha,'Tag','bsFitSelPoints');
for i=1:numel(hBsSelCirc)
    [~,posOfBsSel(i)] = min( abs(hl_splFit.XData - hBsSelCirc(i).XData) );
end

hPeakSelCirc = findall(ha,'Tag','peakFitSelPoints');
[~,posOfPeakSel] = min( abs(hl_splFit.XData - hPeakSelCirc.XData) );

hExpSelCirc = findall(ha,'Tag','expFitSelPoints');
for i=1:numel(hExpSelCirc)
    [~,posOfExpFitSel(i)] = min( abs(hl_splFit.XData - hExpSelCirc(i).XData) );
end

selectedEvent = getappdata(analysisFig,'selectedEvent');
selectedEvent.bsFitS = min(posOfBsSel);
selectedEvent.bsFitE = max(posOfBsSel);

selectedEvent.peakFitPos = posOfPeakSel;
selectedEvent.peakFitVal = hl_splFit.YData(posOfPeakSel);

selectedEvent.expFitS = min(posOfExpFitSel);
selectedEvent.expFitE = max(posOfExpFitSel);
setappdata(analysisFig,'selectedEvent',selectedEvent);

% do fitting and analyzing of profile
findAndAnalyzePeaks([],[],analysisFig)

end

end

