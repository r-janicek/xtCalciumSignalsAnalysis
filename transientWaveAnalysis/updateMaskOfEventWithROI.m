function updateMaskOfEventWithROI(hRect,E)

analysisFig = hRect.Parent.Parent;
selectedEvent = getappdata(analysisFig,'selectedEvent');
h_tempMaskImg = selectedEvent.detectedEvent.h_tempMaskImg;

switch selectedEvent.type
    
    case {'wave','caffeine'}
        tempMask = selectedEvent.detectedEvent.tempMask;
        
    otherwise
        return
end
        
% get mask of rectangle
rectMask = hRect.createMask(h_tempMaskImg);

% update temp mask 
if all(hRect.Color == [ 1 0 0])
    tempMask = tempMask & (~rectMask);
 
elseif all(hRect.Color == [ 0 1 0])
    tempMask = tempMask | rectMask;
    
end
   
% update image
h_tempMaskImg.CData = tempMask;

% save
selectedEvent.detectedEvent.tempMask = tempMask;
setappdata(analysisFig,'selectedEvent',selectedEvent)


end

