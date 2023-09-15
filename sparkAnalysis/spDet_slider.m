function spDet_slider(h_Obj,~,mainFig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

val = round(get(h_Obj,'Value'));
txt = sprintf('watershed sensitivity %d%%',val);
detectEventsPB = getappdata(mainFig,'lastPressedPusbutton'); 
hObjs = getappdata(mainFig,'hObjs');


switch detectEventsPB.Tag
    
    case 'findEvents'
        
        % eventsDetection(detectEventsPB,[],mainFig)
        
    otherwise
        
        profileAnalysis = getappdata(mainFig,'profileAnalysis');
        h_rect_prof_pos = getPosition(profileAnalysis.h_rect_prof);
        plotROIprofile(h_rect_prof_pos,mainFig)
               
end

set(hObjs.txt_spDet,'String',txt)

end

