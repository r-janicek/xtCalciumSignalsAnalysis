function axesTransChFunction(hO, E)
switch E.Button
    case 1 % create window
        % check if window exist
        hf = findall(0, 'Tag','2CH');
        mainFig = findall(0, 'Tag','mainFig');
        hObjs = getappdata(mainFig,'hObjs');
        imgData = getappdata(mainFig,'imgData');
        if isempty(hf)          
            hf = figure('Name','zoom trans channel', ...
                'Tag','2CH', 'Units','normalized', ...
                'Outerposition',[mainFig.Position(1) ...
                                 0.5 ...
                                 mainFig.Position(3) ...
                                 0.4]);
            ax = axes('Parent',hf);
            set(ax,'Position',[hObjs.ax_img.Position(1) ...
                               0.1 ...
                               hObjs.ax_img.Position(3) ...
                               0.89])
            % show image
            % apply crop
            if isfield(imgData, 'cropROIpos')
                if ~isempty(imgData.cropROIpos)
                    r_s = max(min(round(imgData.cropROIpos(2)), size(imgData.imgDataXTtrans,1)), 1);
                    r_e = max(min(round(imgData.cropROIpos(2)+imgData.cropROIpos(4)), size(imgData.imgDataXTtrans,1)), 1);
                    c_s = max(min(round(imgData.cropROIpos(1)/imgData.pxSzT), size(imgData.imgDataXTtrans,2)), 1);
                    c_e = max(min(round((imgData.cropROIpos(1)+imgData.cropROIpos(3))/imgData.pxSzT), size(imgData.imgDataXTtrans,2)), 1);
                else
                    r_s = 1;
                    r_e = size(imgData.imgDataXTtrans,1);
                    c_s = 1;
                    c_e = size(imgData.imgDataXTtrans,2);
                end
            else
                r_s = 1;
                r_e = size(imgData.imgDataXTtrans,1);
                c_s = 1;
                c_e = size(imgData.imgDataXTtrans,2);
            end
            croppedXTtransImg = imgData.imgDataXTtrans(r_s:r_e,c_s:c_e);
            image(imgData.imgDataXTtrans, ...
                'XData',[0 max(imgData.t)], ...
                'YData',[1 size(imgData.imgDataXTtrans,1)], ...
                'CDataMapping','scaled', 'Parent',ax);
            % set(ax,'XTick',[], 'YTick',[], 'PickableParts','all')
            set(ax, 'PickableParts','all')
            xlabel(ax, 'time (ms)')
            % set axes limits so it showing the same as main image axes 
            % (in case slider is activated)
            set(ax, 'XLim',hObjs.ax_img.XLim, ...
                'YLim',hObjs.ax_img.YLim)
        else
            % Figure exists so bring Figure to the focus
            figure(hf);
        end
        
        
    case 3 %delete window
        delete(findobj(0, 'Type','figure', ...
            'Tag','2CH', 'Name','zoom trans channel'))
        
end

end

