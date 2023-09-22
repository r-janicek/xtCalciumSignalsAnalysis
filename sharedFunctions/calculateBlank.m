function calculateBlank(st, ~, mainFig, h_rect_blank)
% calculate blank value from xt image; or by setting value in edit box
style = get(st,'Style');
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');

switch style
    case 'pushbutton'
        % get mask of image area to calculate blank
        bw_blank_crop = createMask(h_rect_blank);
        % delete blank ROI
        delete(h_rect_blank)
        % calculate blank
        blank = mean(imgData.imgDataXTfluoF(bw_blank_crop), 'all');
                      
    case 'edit'
        blank = str2double(get(st,'String'));       
end
  
set(hObjs.h_edit_Blank, 'String',num2str(round(blank,2)))
set(hObjs.h_push_BlankROI, 'String','set ROI for blank', ...
    'FontWeight','normal')
% save data
imgData.blank = blank;
setappdata(mainFig,'imgData',imgData);

end

