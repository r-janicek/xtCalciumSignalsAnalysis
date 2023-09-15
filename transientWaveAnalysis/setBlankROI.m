function setBlankROI(~,~,mainFig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

pxSzT = imgData.pxSzT;
ax_img = hObjs.ax_img;
imgDataXTfluoF = imgData.imgDataXTfluoF;

% blank ROI
h = 10;
if h>=size(imgDataXTfluoF,2)
    h = size(imgDataXTfluoF,2)-1;
end
h_rect_blank = imrect(ax_img,[0 1 size(imgDataXTfluoF,2)*pxSzT h]);
setColor(h_rect_blank,'g')

% constrains for ROI
fcn = makeConstrainToRectFcn('imrect',get(ax_img,'XLim'),get(ax_img,'YLim'));
setPositionConstraintFcn(h_rect_blank,fcn); 
setResizable(h_rect_blank,1);

set(hObjs.h_push_BlankROI,'String','get blank')
set(hObjs.h_push_BlankROI,'Callback',{@calculateBlank,mainFig,h_rect_blank})
set(hObjs.h_push_BlankROI,'FontWeight','bold')

end

