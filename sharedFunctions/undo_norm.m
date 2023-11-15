function undo_norm(~,~,mainFig)
% start over with normalization process

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
    % delete patches texts
    delete(findobj(ax_img,'Type','Text'))
    
    set(whole_img_h, 'CData',imgData.imgDataXTfluoFN)   
    set(get(ax_img, 'Ylabel'),'String','x (pixels) [filtered]')   
    clim(ax_img, [imgData.z_min_img imgData.z_max_img])
    hObjs.ax_img.Toolbar.Visible = 'off';

    if isfield(hObjs, 'ax_img_sparks')
        cla(hObjs.ax_img_sparks)
        hObjs.ax_img_sparks.Toolbar.Visible = 'off';
    end
    cla(hObjs.ax_prof)
    
    set(mainFig,'WindowKeyPressFcn',[]);
    set(mainFig,'WindowButtonUpFcn',[]);
    set(mainFig,'WindowButtonMotionFcn',[]);
     
    % restar table of detected profiles
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
                delete(findall(hObjs.ax_img, 'Type', 'rectangle'))
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

    % show not normalized profile of image
    plot(imgData.t, mean(imgData.imgDataXTfluoFN,1), ...
        'Parent',hObjs.ax_prof, 'Tag','wholeCellProfile')
    set(hObjs.ax_prof, 'FontSize',14, ...
        'Xlim',[imgData.t(1) imgData.t(end)],...
        'YLim', getAxisLimits(mean(imgData.imgDataXTfluoFN,1), 5))
    set(get(hObjs.ax_prof,'Xlabel'), 'String','t (ms)', ...
        'FontWeight','bold')
    set(get(hObjs.ax_prof,'Ylabel'), 'String','Fluorescence (F)', ...
        'FontWeight','bold')
    hObjs.ax_prof.Toolbar.Visible = 'off';
    % disable built-in axes interactions
    arrayfun(@(x) disableDefaultInteractivity(x), ...
        findall(mainFig, 'Type', 'Axes'))
    
    % set up main windows
    switch getappdata(mainFig,'analysisType')
        case {'spark recovery ryanodine', 'spark detection'}
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
            
        case 'transients & waves'
            % reset buttons for ROIs and table
            delete(findobj(ax_img,'Tag','imrect'))
            set(hObjs.h_push_eventImgROI, ...
                'String','<html> <p align="center"> ROI to select image <br> of event <html>')
            set(hObjs.h_push_eventImgROI, 'FontWeight','normal', ...
                'Callback',{@selectEventImage, mainFig})
            set(hObjs.h_table_eventsImgs,'Data',repmat({'----'},[5,1]))

            % set up window
            set(hObjs.h_pb_norm,'Enable','on')
            set(hObjs.h_pb_norm,'String','normalize data')
            set(hObjs.h_pb_norm,'Callback',{@norm_data, mainFig})
            set(hObjs.h_pb_norm,'FontWeight','normal')

            set(hObjs.check_doBsFit,'Enable','on','Value',0)

            set(hObjs.h_pb_GetBlank,'Enable','on')
            set(hObjs.h_push_BlankROI,'Enable','on')
            set(hObjs.h_edit_Blank,'Enable','on')
            %set(hObjs.h_edit_Blank,'String','nan')

            set(hObjs.h_pb_crop,'Enable','on')

            % mouse baseline selection tool
            set(hObjs.ax_prof, 'buttondownfcn',@mouseSetMaskFcn)
            % remove analysis
            if isfield(getappdata(mainFig),'selectedROIs')
                rmappdata(mainFig,'selectedROIs')
            end
    end
    
else
    return
    
end

end

