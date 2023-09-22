function spDet_slider(h_Obj, ~, mainFig)
% slider controlling watershed sensitivity when detecting sparks or

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
% change text of slider
set(hObjs.txt_spDet,'String',txt)

end

