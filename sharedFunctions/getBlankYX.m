function getBlankYX(~, ~, mainFig, h_rect, hf)
% calculate blank value from xy image of whole cell
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
% get data of selected xy image of cell
h_imgXY = findobj(hf,'Type','image');
xyImg = h_imgXY.CData;
% get mask of image area to calculate blank
bw_blank_crop = createMask(h_rect);
% delete blank ROI
delete(h_rect)
% calculate blank
blank = mean(xyImg(bw_blank_crop), 'all');
% save data
imgData.blank = blank;
set(hObjs.h_edit_Blank, 'String',num2str(round(blank,2)))
setappdata(mainFig,'imgData',imgData);
% close figure
close(hf)

end

