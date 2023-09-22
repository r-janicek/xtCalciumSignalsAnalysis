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
    % delete test of patches
    delete(findobj(ax_img,'Type','Text'))
    
    set(whole_img_h,'CData',imgData.imgDataXTfluoFN)   
    set(get(ax_img,'Ylabel'),'String','x (pixels) [filtered]')   
    caxis(ax_img,[imgData.z_min_img imgData.z_max_img])
    
    cla(hObjs.ax_prof)
    
    set(mainFig,'WindowKeyPressFcn',[]);
    set(mainFig,'WindowButtonUpFcn',[]);
    set(mainFig,'WindowButtonMotionFcn',[]);
     
    
    % restar table for detected profiles
    % set(hObjs.h_table_profs,'Data',zeros(5,2))
        
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
    
    % remove last pressed pushbutton
    if isfield(getappdata(mainFig),'lastPressedPusbutton')
        rmappdata(mainFig,'lastPressedPusbutton')
    end
                
    % save data
    imgData.norm_flag = 0;   
    setappdata(mainFig,'imgData',imgData)
    
    % reset buttons for ROIs and table
    delete(findobj(ax_img,'Tag','imrect'))
    set(hObjs.h_push_eventImgROI,'String','<html> <p align="center"> ROI to select image <br> of event <html>')
    set(hObjs.h_push_eventImgROI,'Callback',{@selectEventImage,mainFig})
    set(hObjs.h_push_eventImgROI,'FontWeight','normal')
    set(hObjs.h_table_eventsImgs,'Data',repmat({'----'},[5,1]))
    
    % set up window
    set(hObjs.h_pb_norm,'Enable','on')
    set(hObjs.h_pb_norm,'String','normalize data')
    set(hObjs.h_pb_norm,'Callback',{@norm_data,mainFig})
    set(hObjs.h_pb_norm,'FontWeight','normal')
    
    set(hObjs.check_doBsFit,'Enable','on','Value',0)
    
    set(hObjs.h_pb_GetBlank,'Enable','on')
    set(hObjs.h_push_BlankROI,'Enable','on')
    set(hObjs.h_edit_Blank,'Enable','on')
    %set(hObjs.h_edit_Blank,'String','nan')
    
    set(hObjs.h_pb_crop,'Enable','on')
    
    % show not normalized profile of image
    ax_prof = hObjs.ax_prof;
    plot(imgData.t,mean(imgData.imgDataXTfluoFN,1),'Parent',ax_prof,'Tag','wholeCellProfile')
    set(ax_prof,'Xlim',[imgData.t(1) imgData.t(end)],...
        'YLim',[min(mean(imgData.imgDataXTfluoFN,1)) max(mean(imgData.imgDataXTfluoFN,1))],...
        'FontSize',14)
    set(get(ax_prof,'Xlabel'),'String','t (ms)','FontWeight','bold')
    set(get(ax_prof,'Ylabel'),'String','Fluorescence (F)','FontWeight','bold')
    
    % mouse baseline selection tool
    set(ax_prof,'buttondownfcn',@mouseSetMaskFcn)
    
else
    return
    
end

end

