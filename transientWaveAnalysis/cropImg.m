function cropImg(h_Obj,~,mainFig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

str = h_Obj.String;
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

pxSzT = imgData.pxSzT;
ax_img = hObjs.ax_img;
ax_prof = hObjs.ax_prof;

imgDataXTfluoF = imgData.imgDataXTfluoF;
imgDataXTfluoFN = imgData.imgDataXTfluoFN;
imgDataXTtrans = imgData.imgDataXTtrans;

imgDataXTfluoR = imgData.imgDataXTfluoR;
imgDataXTfluoRN = imgData.imgDataXTfluoRN;

imgDataXT = imgData.imgDataXT;


switch str
    
    case '<html> <p align="center"> set ROI <br> for cropping <html>'
      
        delete(findobj(ax_img,'Type','hggroup','Tag','imrect'))
        h_rect_crop = imrect(ax_img,[size(imgDataXTfluoFN,2)*pxSzT/4 size(imgDataXTfluoFN,1)*0.1 size(imgDataXTfluoFN,2)*pxSzT/2 size(imgDataXTfluoFN,1)*0.8]);
        setColor(h_rect_crop,'r')
        set(h_rect_crop,'DisplayName','croppImg')
        % constrains for ROI
        fcn = makeConstrainToRectFcn('imrect',get(ax_img,'XLim'),get(ax_img,'YLim'));
        setPositionConstraintFcn(h_rect_crop,fcn);
        setResizable(h_rect_crop,1);
        
        set(hObjs.h_pb_crop,'String','crop')
        set(hObjs.h_pb_crop,'FontWeight','bold')
        
        setappdata(mainFig,'h_rect_crop',h_rect_crop);
                    
    case 'crop'
                
        h_rect_crop = getappdata(mainFig,'h_rect_crop');
        
        ROI_pos = getPosition(h_rect_crop);
        pos = [ROI_pos(1)/pxSzT ROI_pos(2) ROI_pos(3)/pxSzT ROI_pos(4)];
        pos(pos<0)=0;
        delete(h_rect_crop)
                
        [r_t(:,1),r_x(:,1)] = rect2ind(pos);
      
        if length(r_t)>size(imgDataXTfluoFN,2)           
            r_t = r_t(1:size(imgDataXTfluoFN,2),1);            
        end
        
        if length(r_x)>size(imgDataXTfluoFN,1)           
            r_x = r_x(1:size(imgDataXTfluoFN,1),1);            
        end
        
        % crop data
        imgDataXTfluoR = imgDataXTfluoR(r_x,r_t);
        imgDataXTfluoRN = imgDataXTfluoRN(r_x,r_t);
        imgDataXTfluoF = imgDataXTfluoF(r_x,r_t);
        imgDataXTfluoFN = imgDataXTfluoFN(r_x,r_t);
        try
            imgDataXTtrans = imgDataXTtrans(r_x,r_t);
        catch
            imgDataXTtrans = [];  
        end
            
        imgDataXT = cellfun(@(x) x(r_x,r_t),imgDataXT,'UniformOutput',0);
        
        % recalculate new t axis
        t = imgData.t;
        t = t(r_t(:,1));
        crop_s_t = t(1); 
        t = t - t(1);
        
        if strcmp(getappdata(mainFig,'analysisType'),'spark recovery photolysis')
            
            % recalculate new positions of photolytic pulses
            [sTPP,eTPP,~,~,~,~] = ...
                loadPhotolysisPositions(imgDataXTfluoF,imgDataXTtrans,...
                t,pxSzT,imgData.pxSzX,...
                imgData.TPPpointPos,imgData.scanLinePos,ax_img,0);
            
        else         
            sTPP = [];
            eTPP = [];                       
        end
        
        % correct positions of lines in cropped image showing photolytic pulses
        lines_TPP = findobj(ax_img,'Type','Line');
        rect_TPP = findobj(ax_img,'Type','Rectangle');
        
        if ~isempty(lines_TPP)
            arrayfun(@(x) set(x,'XData',x.XData - crop_s_t), lines_TPP, 'UniformOutput', false)
            arrayfun(@(x) set(x,'YData',[1 size(imgDataXTfluoFN,1)]), lines_TPP, 'UniformOutput', false)
        end
        
        if ~isempty(rect_TPP)   
            
            arrayfun(@(x) set(x,'Position',x.Position - [crop_s_t  ROI_pos(2) 0 0]), rect_TPP, 'UniformOutput', false)
                                   
        end
        
        set(findobj(ax_img,'Type','Image'),'CData',imgDataXTfluoFN,'XData',[0 max(t)],'YData',[1 size(imgDataXTfluoFN,1)])
        set(ax_img,'XLim',[0 max(t)],'YLim',[1 size(imgDataXTfluoFN,1)])
        
        set(ax_prof.Children,'XData',t,'YData',mean(imgDataXTfluoFN,1))
        set(ax_prof,'YLim',[min(mean(imgDataXTfluoFN,1)) max(mean(imgDataXTfluoFN,1))])
        
        % save data
        imgData.imgDataXT = imgDataXT;
        imgData.imgDataXTfluoR = imgDataXTfluoR;
        imgData.imgDataXTfluoRN = imgDataXTfluoRN;
        imgData.imgDataXTfluoF = imgDataXTfluoF;
        imgData.imgDataXTfluoFN = imgDataXTfluoFN;
        imgData.imgDataXTtrans = imgDataXTtrans;
        imgData.t = t;
        imgData.crop_s_t = crop_s_t;
        imgData.s_TPP = sTPP;
        imgData.e_TPP = eTPP;
        
        
        % save crop roi position
        if isfield(imgData,'cropROIpos')
            cropROIpos = imgData.cropROIpos;   
            imgData.cropROIpos = [cropROIpos(1)+ROI_pos(1) cropROIpos(2)+ROI_pos(2) ROI_pos(3) ROI_pos(4)];
        else
            imgData.cropROIpos = ROI_pos;
        end
                       
        setappdata(mainFig,'imgData',imgData);
        
        % set up pushbutton for cropping
        set(hObjs.h_pb_crop,'String','<html> <p align="center"> set ROI <br> for cropping <html>')
        set(hObjs.h_pb_crop,'FontWeight','normal')
               
end

end

