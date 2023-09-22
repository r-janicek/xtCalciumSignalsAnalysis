function norm_data(~,~,mainFig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
keyboard
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');

pxSzT = imgData.pxSzT;
blank = str2double(get(hObjs.h_edit_Blank,'String'));
if isnan(blank)
    
    hObjs.h_edit_Blank.String = '0';
    blank = str2double(get(hObjs.h_edit_Blank,'String'));
%     warndlg('get blank');
%     return
end

ax_img = hObjs.ax_img;
imgDataXTfluoF = imgData.imgDataXTfluoF;
ax_img_yLim = get(ax_img,'YLim');

% F0 ROI
h_rect_F0 = imrect(ax_img,[0 ax_img_yLim(1) size(imgDataXTfluoF,2)*pxSzT/8 ax_img_yLim(2)]);
setColor(h_rect_F0,'g')

% constrains for ROI
fcn = makeConstrainToRectFcn('imrect',get(ax_img,'XLim'),get(ax_img,'YLim'));
setPositionConstraintFcn(h_rect_F0,fcn); 
setResizable(h_rect_F0,1);

set(hObjs.h_pb_norm,'String','normalize (F/F0)')
set(hObjs.h_pb_norm,'Callback',{@norm_calc,mainFig,h_rect_F0})
set(hObjs.h_pb_norm,'FontWeight','bold')


 