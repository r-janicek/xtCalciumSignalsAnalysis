function cropImg(h_Obj,~,mainFig)

str = h_Obj.String;
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;
ax_img = hObjs.ax_img;
ax_prof = hObjs.ax_prof;

imgDataXTfluoF = imgData.imgDataXTfluoF;
imgDataXTfluoFN = imgData.imgDataXTfluoFN;
imgDataXTtrans = imgData.imgDataXTtrans;

imgDataXTfluoR = imgData.imgDataXTfluoR;
imgDataXTfluoRN = imgData.imgDataXTfluoRN;

imgDataXT = imgData.imgDataXT;

switch str
    case '<html> <p align="center"> set ROI to crop <html>'
        
        delete(findobj(ax_img,'Type','hggroup','Tag','imrect'))
        % suggest size of cropping rectangle
        % get 10th percentile from xy image - estimation of background
        if ~isempty(imgData.imgDataXY)
            
            bg_est = double(prctile(imgData.imgDataXY{1}(:),10));
            if prctile(imgData.wholeImgFluoXT(:),1) < bg_est*1.10
                % get image data
                mImg = imgData.imgDataXTfluoFN - bg_est*1.10;
                mImg(mImg>0) = 1;
                mImg(mImg<=0) = 0;
              
                try
                    % find biggest rectangle, which is in image area
                    img_CC_event = bwconncomp(logical(mImg),8);
                    
                    %calculate properties of all CC (regions)
                    statEvents = regionprops(img_CC_event,imgData.imgDataXTfluoFN, ...
                        'SubarrayIdx', 'Area','BoundingBox');
                    [~,indMaxArea] = max([statEvents.Area]);
                    h_rect_crop_newPos = statEvents(indMaxArea).BoundingBox;
                    h_rect_crop_newPos(1) = h_rect_crop_newPos(1)*imgData.pxSzT;
                    h_rect_crop_newPos(3) = h_rect_crop_newPos(3)*imgData.pxSzT - ...
                        h_rect_crop_newPos(3)*imgData.pxSzT*0.01;
                    
                    h_rect_crop = drawrectangle(ax_img, ...
                        'Position',h_rect_crop_newPos);
                    
                catch
                    h_rect_crop = drawrectangle(ax_img, ...
                        'Position',...
                        [ax_img.XLim(1)+2*imgData.pxSzT ax_img.YLim(1)+2.5 ...
                        ax_img.XLim(2)-2*imgData.pxSzT ax_img.YLim(2)-2.5]);
                end
            else
                xr = size(imgDataXTfluoFN,2)*pxSzT/4;
                yr = size(imgDataXTfluoFN,1)*0.1;
                wr = size(imgDataXTfluoFN,2)*pxSzT/2;
                hr = size(imgDataXTfluoFN,1)*0.8;
                h_rect_crop = drawrectangle(ax_img, ...
                    'Position',[xr yr wr hr]);
            end
            
        else
            xr = size(imgDataXTfluoFN,2)*pxSzT/4;
            yr = size(imgDataXTfluoFN,1)*0.1;
            wr = size(imgDataXTfluoFN,2)*pxSzT/2;
            hr = size(imgDataXTfluoFN,1)*0.8;
            h_rect_crop = drawrectangle(ax_img, ...
                'Position',[xr yr wr hr]);
        end
        % setup roi
        % try to set postition of crop ROI 
        if imgData.t(end) > 30000
            h_rect_crop.Position = [
                6500 ...
                size(imgDataXTfluoFN,1)*0.1 ...
                30000 ...
                size(imgDataXTfluoFN,1)*0.8];
        end
        h_rect_crop.Color = 'r';
        h_rect_crop.Tag = 'croppImg';

        set(hObjs.h_pb_crop,'String','crop')
        set(hObjs.h_pb_crop,'FontWeight','bold')
        
        setappdata(mainFig,'h_rect_crop',h_rect_crop);
        
    case 'crop'
        % get position of crop ROI        
        h_rect_crop = getappdata(mainFig,'h_rect_crop');
        bw_crop = createMask(h_rect_crop);
        r_t = find(any(bw_crop, 1));
        r_x = find(any(bw_crop, 2));
        ROI_pos = h_rect_crop.Position;
        % delete crop ROI
        delete(h_rect_crop)
                
        if length(r_t)>size(imgDataXTfluoFN, 2)           
            r_t = r_t(1:size(imgDataXTfluoFN,2),1);            
        end
        if length(r_x)>size(imgDataXTfluoFN, 1)           
            r_x = r_x(1:size(imgDataXTfluoFN, 1), 1);            
        end
        
        % crop data
        imgDataXTfluoR = imgDataXTfluoR(r_x,r_t);
        imgDataXTfluoRN = imgDataXTfluoRN(r_x,r_t);
        imgDataXTfluoF = imgDataXTfluoF(r_x,r_t);
        imgDataXTfluoFN = imgDataXTfluoFN(r_x,r_t);
        % recalculate new t axis
        t = imgData.t;
        t = t(r_t);
        crop_s_t = t(1); 
        t = t - t(1);
        % crop trans data
        try
            imgDataXTtrans = imgDataXTtrans(r_x,r_t);
            set(findobj(hObjs.h_ax_transCh, 'Type','Image'), ...
                'CData',imgDataXTtrans, ...
                'XData',[0 max(t)], ...
                'YData',[1 size(imgDataXTtrans,1)])
            set(hObjs.h_ax_transCh, 'XLim',[0 max(t)], ...
                'YLim',[1 size(imgDataXTtrans,1)])
        catch
            imgDataXTtrans = [];  
        end
        imgDataXT = cellfun(@(x) x(r_x,r_t), imgDataXT, ...
            'UniformOutput',0);
        
        % set image axes limits
        set(findobj(ax_img, 'Type','Image'), ...
            'CData',imgDataXTfluoFN, ...
            'XData',[0 max(t)], ...
            'YData',[1 size(imgDataXTfluoFN,1)])
        set(ax_img, 'XLim',[0 max(t)], ...
            'YLim',[1 size(imgDataXTfluoFN,1)])
        % set profile
        % check if there is baseline mask or fit
        if isfield(imgData, 'baselineM')
            % delete previous baseline fit line
            delete(findobj(hObjs.ax_prof, 'Type','Line', ...
                'Tag','baselineFit'))
            % delete previous baseline mask line
            delete(findobj(hObjs.ax_prof, 'Type','Line', ...
                '-regexp','Tag','Mask'))
            % setup length of baseline mask
            imgData.baselineM = imgData.baselineM(r_t);
        end
        set(ax_prof.Children, 'XData',t, 'YData',mean(imgDataXTfluoFN,1), ...
            'Tag','wholeCellProfile')
        set(ax_prof, ...
            'YLim',getAxisLimits(mean(imgDataXTfluoFN,1), 5))
        % try to clear img sparks axes
        if isfield(hObjs, 'ax_img_sparks')
            cla(hObjs.ax_img_sparks); 
        end
        % save data
        imgData.imgDataXT = imgDataXT;
        imgData.imgDataXTfluoR = imgDataXTfluoR;
        imgData.imgDataXTfluoRN = imgDataXTfluoRN;
        imgData.imgDataXTfluoF = imgDataXTfluoF;
        imgData.imgDataXTfluoFN = imgDataXTfluoFN;
        imgData.imgDataXTtrans = imgDataXTtrans;
        imgData.t = t;
        imgData.crop_s_t = crop_s_t;
       
        % save crop roi position
        if isfield(imgData,'cropROIpos')
            cropROIpos = imgData.cropROIpos;   
            imgData.cropROIpos = [cropROIpos(1)+ROI_pos(1), ... 
                                  cropROIpos(2)+ROI_pos(2), ...
                                  ROI_pos(3), ... 
                                  ROI_pos(4)];
        else
            imgData.cropROIpos = ROI_pos;
        end
        % save data to gui               
        setappdata(mainFig,'imgData',imgData);
        
        % set up pushbutton for cropping
        set(hObjs.h_pb_crop, ...
            'String','<html> <p align="center"> set ROI to crop <html>')
        set(hObjs.h_pb_crop,'FontWeight','normal')
        
        % set width of profile ROI
        h_ROI_prof = hObjs.ax_img.YLim(2)-hObjs.ax_img.YLim(1)-4;
        if mod(h_ROI_prof,2)==0
            h_ROI_prof = h_ROI_prof+1;
        end
 
        try 
            set(hObjs.h_edit_ROI,'String',num2str(h_ROI_prof))
        catch             
        end

        % show/calculate baseline mask
        if strcmp(mainFig.Name, 'transients and waves analysis')
            checkboxFcn(hObjs.check_doBsFit, [], mainFig)
        end
end

end

