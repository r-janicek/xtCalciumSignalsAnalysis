function scaleLineChange(listenerProp, mainFig)
% change scale string
imgData = getappdata(mainFig, 'imgData');
hObjs = getappdata(mainFig, 'hObjs');
% update scale bar text
sc_num_img = (diff(hObjs.ax_img.YLim) * imgData.pxSzX)/2; % in um
hObjs.h_txt_scale_img.String = [sprintf('%0.2f',sc_num_img),' \mum'];
try
    sc_num_sparks = (diff(hObjs.ax_img_sparks.YLim) * imgData.pxSzX)/2; % in um
    hObjs.h_txt_scale_sparks.String = [sprintf('%0.2f',sc_num_sparks),' \mum'];
catch
end
end

