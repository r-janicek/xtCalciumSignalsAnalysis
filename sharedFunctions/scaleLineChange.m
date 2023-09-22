function scaleLineChange(listenerProp,mainFig)

imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

% set and add interactive scale
lineH = hObjs.ax_img.Position(4)/hObjs.ax_sc.Position(4); % line height as ratio between axes for image and invisible axes for scale line

sc_num = (diff(hObjs.ax_img.YLim)*imgData.pxSzX)/lineH; % in um
hObjs.h_txt_scale.String = [sprintf('%0.2f',sc_num),' \mum'];


end

