function s = keyPressProfileROImovement(mainFig,s,h_rect,sz_img)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

try
    p = getPosition(h_rect);
    
    pxSize_t = getappdata(mainFig,'pxSize_t');
    
    h = sz_img(1);
    w = sz_img(2)*pxSize_t;
    
    switch s.Key
        
        case 'rightarrow'
            if p(1)+p(3) < w
                setPosition(h_rect,[p(1)+pxSize_t p(2) p(3) p(4)])
            end
            
        case 'leftarrow'
            if p(1) > 0
                setPosition(h_rect,[p(1)-pxSize_t p(2) p(3) p(4)])
            end
            
        case 'uparrow'
            
            if p(2) > 0
                setPosition(h_rect,[p(1) p(2)-1 p(3) p(4)])
            end
        case 'downarrow'
            if p(2)+p(4) < h
                setPosition(h_rect,[p(1) p(2)+1 p(3) p(4)])
            end
        
        otherwise
            return
                        
    end
       
    plotROIprofile(getPosition(h_rect),mainFig)
    
catch
    
end

end

