function checkRectVisibility(h_obj,~,mainFig)

val = get(h_obj,'Value');
sparkDetection = getappdata(mainFig,'sparkDetection');
hObjs = getappdata(mainFig,'hObjs');
rectTxt = findall(hObjs.ax_img,'Tag','detectedEventRecText');

if isfield(sparkDetection,'detectedEventsRec')
    
    if val == 1
        set(sparkDetection.detectedEventsRec,'Visible','on')
        set(sparkDetection.detectedEventsMask,'Visible','on')
        set(rectTxt,'Visible','on')
    else
        set(sparkDetection.detectedEventsRec,'Visible','off')
        set(sparkDetection.detectedEventsMask,'Visible','off')
        set(rectTxt,'Visible','off')
    end
    
else
    set(h_obj,'Value',0)
    return   
end

end

