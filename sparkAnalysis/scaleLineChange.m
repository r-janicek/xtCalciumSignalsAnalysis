function scaleLineChange(listenerProp,mainFig)

imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

% set and add interactive scale
sc_num = (diff(hObjs.ax_img_sparks.YLim)*imgData.pxSzX)/2; % in um
hObjs.h_txt_scale.String = [sprintf('%0.2f',sc_num),' \mum'];


end

