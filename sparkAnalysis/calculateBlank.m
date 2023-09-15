function calculateBlank(st,~,mainFig,h_rect_blank)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

style = get(st,'Style');
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');

pxSzT = imgData.pxSzT;
imgDataXTfluoF = imgData.imgDataXTfluoF;

switch style
    
    case 'pushbutton'
        
        %get pos of blank ROI
        ROI_pos = getPosition(h_rect_blank);
        pos = [ROI_pos(1)/pxSzT ROI_pos(2) ROI_pos(3)/pxSzT ROI_pos(4)];
        pos(pos<0)=0;
         
        delete(h_rect_blank)
        
        %calculate blank
        
        [r_t(:,1),r_x(:,1)] = rect2ind(pos);
      
        if length(r_t)>size(imgDataXTfluoF,2)           
            r_t = r_t(1:size(imgDataXTfluoF,2),1);            
        end
        
        if length(r_x)>size(imgDataXTfluoF,1)           
            r_x = r_x(1:size(imgDataXTfluoF,1),1);            
        end
             
        crop = imgDataXTfluoF(r_x,r_t);       
        blank = mean(crop(:));
                      
    case 'edit'
        blank = str2double(get(st,'String'));
        
end

imgData.blank = blank;

set(hObjs.h_edit_Blank,'String',num2str(round(blank,2)))
set(hObjs.h_push_BlankROI,'String','set ROI for blank','FontWeight','normal')

setappdata(mainFig,'imgData',imgData);

end

