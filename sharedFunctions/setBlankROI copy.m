function setBlankROI(~,~,mainFig)

imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

pxSzT = imgData.pxSzT;
ax_img = hObjs.ax_img;
imgDataXTfluoF = imgData.imgDataXTfluoF;

% blank ROI
h = 10;
if h >= size(imgDataXTfluoF,2)
    h = size(imgDataXTfluoF,2)-1;
end
h_rect_blank = drawrectangle(ax_img, ...
    'Position',[0 1 size(imgDataXTfluoF,2)*pxSzT h], ...
    'Color','g');
% setup pushbutton
set(hObjs.h_push_BlankROI, 'String','get blank', 'FontWeight','bold', ...
    'Callback',{@calculateBlank,mainFig,h_rect_blank})
end

