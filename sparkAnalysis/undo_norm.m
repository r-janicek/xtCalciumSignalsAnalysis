function undo_norm(~,~,mainFig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs'); 

if imgData.norm_flag==1

    % set up data
    imgData.imgDataXTfluoRN = imgData.imgDataXTfluoR;
    imgData.imgDataXTfluoFN = imgData.imgDataXTfluoF;
        
    imgData.z_max_img = max(imgData.imgDataXTfluoF(:));
    imgData.z_min_img = min(imgData.imgDataXTfluoF(:));
     
    ax_img = hObjs.ax_img;
    whole_img_h = findobj(ax_img,'Type','Image');
    % delete patches of profiles
    delete(findobj(ax_img,'Type','Patch'))
       
    set(whole_img_h,'CData',imgData.imgDataXTfluoFN)   
    set(get(ax_img,'Ylabel'),'String','x (pixels) [filtered]')   
    caxis(ax_img,[imgData.z_min_img imgData.z_max_img])
    
    cla(hObjs.ax_img_sparks)
    cla(hObjs.ax_prof)
    
    set(mainFig,'WindowKeyPressFcn',[]);
    set(mainFig,'WindowButtonUpFcn',[]);
    set(mainFig,'WindowButtonMotionFcn',[]);
     
    
    % restar table for detected profiles
    if isfield(getappdata(mainFig),'profileAnalysis')
        set(hObjs.h_table_profs,'Data',zeros(5,2))
    end
    
    % delete imrect of roi for detecting profiles if exist,
    % also remove analysis
    if isfield(getappdata(mainFig),'profileAnalysis')
        profileAnalysis = getappdata(mainFig,'profileAnalysis');
        try
            delete(profileAnalysis.h_rect_prof)
        catch
        end
        
        rmappdata(mainFig,'profileAnalysis')
    end

    % remove previous spark detection analysis
    if isfield(getappdata(mainFig),'sparkDetection')
        sparkDetection = getappdata(mainFig,'sparkDetection');
        try
            if isfield(sparkDetection,'detectedEventsRec')
                delete(findall(hObjs.ax_img,'Tag','detectedEventRecText'))
                delete(sparkDetection.detectedEventsRec)
            end
        catch
        end
        rmappdata(mainFig,'sparkDetection')
    end
    
    % remove last pressed pushbutton
    if isfield(getappdata(mainFig),'lastPressedPusbutton')
        rmappdata(mainFig,'lastPressedPusbutton')
    end
                
    % save data
    imgData.norm_flag = 0;   
    setappdata(mainFig,'imgData',imgData)
    
    % set up window
    switch getappdata(mainFig,'analysisType')
        case 'spark detection'
            set(hObjs.h_pb_norm,'Enable','off')
        otherwise
            set(hObjs.h_pb_norm,'Enable','on')
    end
    set(hObjs.h_pb_norm,'String','normalize data')
    set(hObjs.h_pb_norm,'Callback',{@norm_data,mainFig})
    set(hObjs.h_pb_norm,'FontWeight','normal')
    
    set(hObjs.h_pb_GetBlank,'Enable','on')
    set(hObjs.h_push_BlankROI,'Enable','on')
    set(hObjs.h_edit_Blank,'Enable','on')
    %set(hObjs.h_edit_Blank,'String','nan')
    
    set(hObjs.h_pb_crop,'Enable','on')
   
else
    return
    
end

end

