function scaleLineChange(listenerProp,mainFig)

imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
% update scale bar text
sc_num = (diff(hObjs.ax_img.YLim) * imgData.pxSzX)/2; % in um
hObjs.h_txt_scale_img.String = [sprintf('%0.2f',sc_num),' \mum'];

end

