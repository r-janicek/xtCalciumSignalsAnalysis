function setPxSzManually(hO, E, mainFig)
% set pixel sizes manually in edit box
% get changed pixel size
pxSz = str2double(hO.String);
imgData = getappdata(mainFig, 'imgData');
hObjs = getappdata(mainFig, 'hObjs');
switch hO.Tag
    case 'pxSzX'
        imgData.pxSzX = pxSz;
        delete(getappdata(mainFig,'editScaleListener_sparks'))
        delete(getappdata(mainFig,'editScaleListener_img'))
        % change scale
        hObjs.h_txt_scale_img.String = ...
            [sprintf('%0.2f',(diff(hObjs.ax_img.YLim)*pxSz)/2),' \mum'];
        editScaleListener_img = addlistener(hObjs.ax_img, ...
            'YLim', 'PostSet', ...
            @(varargin)scaleLineChange(varargin,mainFig));
        setappdata(mainFig,'editScaleListener_img',editScaleListener_img)
        if isfield(hObjs, 'ax_img_sparks')
            hObjs.h_txt_scale_sparks.String = ...
                [sprintf('%0.2f',(diff(hObjs.ax_img_sparks.YLim)*pxSz)/2),' \mum'];
            editScaleListener_sparks = addlistener(hObjs.ax_img_sparks, ...
                'YLim', 'PostSet', ...
                @(varargin)scaleLineChange(varargin,mainFig));
            setappdata(mainFig,'editScaleListener_sparks',editScaleListener_sparks)
        end
        
    case 'pxSzT'
        imgData.pxSzT = pxSz;
end
% save changed image data
setappdata(mainFig, 'imgData', imgData)
end