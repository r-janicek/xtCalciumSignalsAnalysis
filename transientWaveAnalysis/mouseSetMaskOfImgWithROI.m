function mouseSetMaskOfImgWithROI(hROI,E,analysisFig)
keyboard
switch E.Button
    case 1   
        hROI.EdgeColor = 'r';
        
    case 3
        hROI.EdgeColor = 'g';
end

try
%     % get image mask
%     selectedEvent = getappdata(analysisFig,'selectedEvent');
%     h_tempMaskImg = selectedEvent.detectedEvent.h_tempMaskImg;
%     
%     
%     cp = hO.Parent.CurrentPoint;
%     
%     xinit = cp(1,1);
%     yinit = cp(1,2);
%     hl = line('XData',xinit,'YData',yinit,...
%         'Marker','none' ,'color','g','LineWidth',5);
    if isempty(analysisFig.WindowButtonMotionFcn)
        analysisFig.WindowButtonMotionFcn = {@wbmcb,hROI};
    else
        analysisFig.WindowButtonMotionFcn = {@wbmcb,hROI};
        analysisFig.WindowButtonUpFcn = {@wbucb,hROI};   
    end
    
catch
    return
end
 

% function executed when move with mouse pointer
function wbmcb(analysisFig,~,hROI)

cp = hROI.Parent.CurrentPoint;

if (cp(1,1) < hROI.Parent.XLim(1)) || (cp(1,1) > hROI.Parent.XLim(2)) || ...
        (cp(1,2) < hROI.Parent.YLim(1)) || (cp(1,2) > hROI.Parent.YLim(2))
    
    hROI.Position = [hROI.Parent.XLim(1), hROI.Parent.YLim(1), ...
        hROI.Position(3), hROI.Position(4)];

else
    hROI.Position = [round(cp(1,1)-hROI.Position(3)/2), ...
        round(cp(1,2)-hROI.Position(4)/2), ...
        hROI.Position(3), hROI.Position(4)];
end

roiPos = hROI.Position;

if ~isempty(analysisFig.WindowButtonMotionFcn)
    
    % centre of ROI
    cx = round(roiPos(1)+roiPos(3)/2);
    cy = round(roiPos(2)+roiPos(4)/2);
    
    % axis of ellipse
    ex = roiPos(3)/2;
    ey = roiPos(4)/2;
    
    % get temp mask
    selectedEvent = getappdata(analysisFig,'selectedEvent');
    tempMask = selectedEvent.detectedEvent.tempMask;
    
    % create ROI mask from ellipse
    [cols_roiMask, rows_roiMask] = meshgrid(1:size(tempMask,2), 1:size(tempMask,1));
    % 2D array
    roiMask = ((rows_roiMask - cy).^2)./ey.^2 ...
        + ((cols_roiMask - cx).^2)./ex.^2 <= 1;
    
    % update temp mask
    if all(hROI.EdgeColor == [ 1 0 0])
        tempMask = tempMask & (~roiMask);
        
    elseif all(hROI.EdgeColor == [ 0 1 0])
        tempMask = tempMask | roiMask;
        
    end
    
    % update image
    selectedEvent.detectedEvent.h_tempMaskImg.CData = tempMask;
    drawnow
    
    % save
    selectedEvent.detectedEvent.tempMask = tempMask;
    setappdata(analysisFig,'selectedEvent',selectedEvent)
end

end

% stop function, when mouse button is released
function wbucb(analysisFig,~,hROI)
% 
 analysisFig.WindowButtonMotionFcn = '';
 analysisFig.WindowButtonUpFcn = '';
% 
% % set mask of events
% hl_x = hl.XData;
% 
% imgData = getappdata(fig,'imgData');
% hObjs = getappdata(fig,'hObjs');
% 
% m = imgData.baselineM;
%  
% xt = imgData.t;
% profLineObj = findall(ha,'Tag','wholeCellProfile');
% prof_t = profLineObj.YData;
% 
% if numel(hl_x)<2
%     return
% end
%     
% % position of start and end of line in points
% [~,s] = min(abs(xt - hl_x(1)));
% [~,e] = min(abs(xt - hl_x(2)));
% 
% if s <= e
%     % mark new area for event fitting 
%     m(s:e,1) = true;
%        
% else
%     % remove mark from selected area
%     m(e:s,1) = false;      
% end
% 
% % delete selection line and marked events line
% delete(hl)
% eventsMaskLine = findobj(ha,'Type','Line','-regexp','Tag','Mask');
% delete(eventsMaskLine)
% 
% % create new line
% try
%     hl_m = line(xt(m),prof_t(m),'Parent',ha,'Color','g',...
%         'LineStyle','none','Marker','.','MarkerSize',20,'LineWidth',1,'Tag','Mask');
%     uistack(hl_m, 'bottom')
% catch
%     keyboard
% end
% 
% % do fitting of baseline
% type = hObjs.popUpMenuBaselineFcn.String{hObjs.popUpMenuBaselineFcn.Value};
% bsFit = fitBaseline(xt,prof_t,type,m,1,ha);
% 
% % save new baseline mask
% imgData.baselineM = m;
% imgData.bsFit = bsFit;
% setappdata(fig,'imgData',imgData);

end



end

