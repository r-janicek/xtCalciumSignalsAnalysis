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
        old_pxSzT = imgData.pxSzT;
        imgData.pxSzT = pxSz;
        % re-create time vector, starting at 0 (ms)
        imgData.t = linspace(0, ...
                             (size(imgData.imgDataXTfluoFN,2)-1)*pxSz, ...
                             size(imgData.imgDataXTfluoFN,2));
        % set up axes
        h_img = findobj(hObjs.ax_img, 'Type','image');
        h_img.XData = [min(imgData.t) max(imgData.t)];
        hObjs.ax_img.XLim = [min(imgData.t) max(imgData.t)];
        h_prof = findobj(hObjs.ax_prof, ...
            {'Type','line'},'-and', ...
            {'-not',{{'-regexp','Tag','Mask'},'-or',{'-regexp','Tag','baseline'}}}); 
        h_prof.XData = imgData.t;
        % set up
        if isfield(imgData, 'crop_s_t')
            imgData.crop_s_t = imgData.crop_s_t*(pxSz/old_pxSzT);
        end

end
% filter data
[imgDataXTfluo, imgData.imgFiltersUsed] = imgFiltering( ...
    imgData.imgDataXTfluoR, imgData.pxSzT, imgData.pxSzX);
% save changed image data
imgData.imgDataXTfluoF = imgDataXTfluo;
imgData.imgDataXTfluoFN =imgDataXTfluo;
setappdata(mainFig, 'imgData', imgData)

end