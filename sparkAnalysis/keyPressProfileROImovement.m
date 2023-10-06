function E = keyPressProfileROImovement(mainFig, E, h_rect, sz_img)
% move profile ROI with keyboard arrows
try
    % get ROI position
    p = h_rect.Position;
    imgData = getappdata(mainFig,'imgData');
    pxSzT = imgData.pxSzT;
    % not to move ROI out of image
    h = sz_img(1);
    w = sz_img(2)*pxSzT;

    switch E.Key
        
        case 'rightarrow'
            if p(1)+p(3) < w
                h_rect.Position = [p(1)+pxSzT p(2) p(3) p(4)];
            end
            
        case 'leftarrow'
            if p(1) > 0
                h_rect.Position = [p(1)-pxSzT p(2) p(3) p(4)];
            end
            
        case 'uparrow'
            if p(2) > 0
                h_rect.Position = [p(1) p(2)-1 p(3) p(4)];
            end
        case 'downarrow'
            if p(2)+p(4) < h
                h_rect.Position = [p(1) p(2)+1 p(3) p(4)];
            end
        
        otherwise
            return             
    end
    % plot profile from moved ROI   
    plotROIprofile(h_rect, [], mainFig)
    
catch
    
end

end

