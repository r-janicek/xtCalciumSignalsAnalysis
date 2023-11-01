function norm_data(~,~,mainFig)
% rectangle to select a area from which mean value of baseline is calulated
% for each time profile
hObjs = getappdata(mainFig, 'hObjs');
imgData = getappdata(mainFig, 'imgData');
pxSzT = imgData.pxSzT;
blank = str2double(get(hObjs.h_edit_Blank, 'String'));
if isnan(blank)
    errordlg('MEASURE BLANK!')
    return
end
ax_img = hObjs.ax_img;
imgDataXTfluoF = imgData.imgDataXTfluoF;
ax_img_yLim = get(ax_img, 'YLim');
% if baseline fitting is selected do not draw rectangle to select area to
% calculate F0
if hObjs.check_doBsFit.Value
    set(hObjs.h_pb_norm, 'String','normalizing')
    set(hObjs.h_pb_norm, 'FontWeight','bold')
    norm_calc([],[],mainFig,[])
else
    % F0 ROI
    h_rect_F0 = drawrectangle(ax_img, 'Color','g',...
        'Position',[0 ...
                    ax_img_yLim(1) ...
                    size(imgDataXTfluoF,2)*pxSzT/8 ...
                    ax_img_yLim(2)]);
    % setup pushbutton
    set(hObjs.h_pb_norm, 'String','normalize')
    set(hObjs.h_pb_norm, 'Callback',{@norm_calc,mainFig,h_rect_F0})
    set(hObjs.h_pb_norm, 'FontWeight','bold')
end


 