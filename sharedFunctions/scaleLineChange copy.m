function scaleLineChange(listenerProp, mainFig)
% change scale string
imgData = getappdata(mainFig, 'imgData');
hObjs = getappdata(mainFig, 'hObjs');
% change text of scale; in um
sc_num = (diff(listenerProp{1,2}.AffectedObject.YLim)*imgData.pxSzX)/2; 
switch listenerProp{1,2}.AffectedObject.Tag
    case 'img'
        hObjs.h_txt_scale_img.String = ...
            [sprintf('%0.2f',sc_num),' \mum'];
    case 'img_sparks'
        hObjs.h_txt_scale_sparks.String = ...
            [sprintf('%0.2f',sc_num),' \mum'];
end

end

