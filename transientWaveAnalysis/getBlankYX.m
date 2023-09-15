function getBlankYX(~, ~,mainFig,h_rect,hf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

% get data of selected xy image of cell
h_imgXY = findobj(hf,'Type','image');
xyImg = h_imgXY.CData;

% get selected area
pos_rect = getPosition(h_rect);
pos_rect(pos_rect<1)=1;

[r_x(:,1),r_y(:,1)] = rect2ind(pos_rect);
crop = xyImg(r_y,r_x);

% calculate blank as mean of selected area
blank = mean(crop(:));

% save data
imgData.blank = blank;

set(hObjs.h_edit_Blank,'String',num2str(round(blank,2)))
setappdata(mainFig,'imgData',imgData);

close(hf)


end

