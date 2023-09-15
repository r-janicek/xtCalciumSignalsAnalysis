function selectEventImage(~,~,mainFig)

hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');
pxSzT = imgData.pxSzT;
ax_img = hObjs.ax_img;
imgDataXTfluoF = imgData.imgDataXTfluoF;
ax_img_yLim = get(ax_img,'YLim');

eventType = hObjs.h_bg_eventSelection.SelectedObject.String;

% ROI
% % delete previous imrect and set wave pushbutton to default state
% delete(findobj(ax_img,'Tag','imrect'))
% set(hObjs.h_push_transROI,'String','<html> <p align="center"> ROI to select image <br> of transients <html>')
% set(hObjs.h_push_transROI,'Callback',{@selectTransientsImage,mainFig})
% set(hObjs.h_push_transROI,'FontWeight','normal')

% create new ROI imrect
wRoi = 2000; % ms
if diff(ax_img.XLim) < wRoi
    wRoi = diff(ax_img.XLim)/2;
end

h_rect = imrect(ax_img,[ax_img.XLim(1) ax_img_yLim(1) wRoi ax_img_yLim(2)]);

switch eventType
    case 'wave'
        setColor(h_rect,'r')
        
    case 'transient'
        setColor(h_rect,'k')
        
    case 'caffeine'
        setColor(h_rect,'g')    
end
    
% constrains for ROI
% use constrains from image, in axes
imgInAxes = findall(ax_img,'Type','image');
try
    fcn = makeConstrainToRectFcn('imrect',imgInAxes.XData,imgInAxes.YData);
catch
    fcn = makeConstrainToRectFcn('imrect',ax_img.XLim,ax_img.YLim);
end
setPositionConstraintFcn(h_rect,fcn); 
setResizable(h_rect,1);

% change pushbutton
set(hObjs.h_push_eventImgROI,'String','<html> <p align="center"> get image <br> of event <html>')
set(hObjs.h_push_eventImgROI,'Callback',{@getImageData,mainFig,h_rect})
set(hObjs.h_push_eventImgROI,'FontWeight','bold')

% disable change of type of event
hObjs.h_b1_eventSelection.Enable = 'off';
hObjs.h_b2_eventSelection.Enable = 'off';
hObjs.h_b3_eventSelection.Enable = 'off';
hObjs.h_b4_eventSelection.Enable = 'off';

end


