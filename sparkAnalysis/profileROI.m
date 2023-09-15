function profileROI(h_obj,~,mainFig)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

str = get(h_obj,'String');
stl = get(h_obj,'Style');

imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
profileAnalysis = getappdata(mainFig,'profileAnalysis');

setappdata(mainFig,'lastPressedPusbutton',h_obj);

imgDataXTfluoR = imgData.imgDataXTfluoR;
ax_prof = hObjs.ax_prof;
ax_img = hObjs.ax_img;

pxSzT = imgData.pxSzT;
t = imgData.t;

switch stl
    
    case 'pushbutton'
                    
        if strcmp('whole-cell profile',str)
                            
            if isfield(profileAnalysis,'h_rect_prof')                
                delete(profileAnalysis.h_rect_prof)                           
            end
            keyboard
            h_rect_prof_pos = [0 1 size(imgDataXTfluoR,2)*pxSzT size(imgDataXTfluoR,1)];
            plotROIprofile(h_rect_prof_pos,mainFig)
            
        else
          
            h = str2double(get(hObjs.h_edit_ROI,'String'));
         
            md = get(ax_img,'YLim');
            
            if isfield(profileAnalysis,'h_rect_prof')                
                delete(profileAnalysis.h_rect_prof)                           
            end
       
            % ROI
            h_rect_prof = imrect(ax_img,[t(1) md(1) t(end) h]);
            setColor(h_rect_prof,'r')
            set(h_rect_prof,'PickableParts','all')
            
            profileAnalysis.h_rect_prof = h_rect_prof;          
            setappdata(mainFig,'profileAnalysis',profileAnalysis);
            
            % constrains for ROI
            fcn = makeConstrainToRectFcn('imrect',get(ax_img,'XLim'),get(ax_img,'YLim'));
            setPositionConstraintFcn(h_rect_prof,fcn);
            setResizable(h_rect_prof,0);
    
            %h_rect_prof_pos = getPosition(h_rect_prof);
            %set(h_rect_prof,'ButtonDownFcn',{@c})
            %h_rect_prof.addNewPositionCallback(@(h_rect_prof_pos) plotProfile(h_rect_prof_pos,main_fig))   
            set(mainFig,'WindowKeyPressFcn',{@keyPressProfileROImovement,h_rect_prof,size(imgDataXTfluoR)});
            set(mainFig,'WindowButtonUpFcn',@mouseProfileRoiMovement);
           
        end
              
    case 'edit'
    
        if isfield(profileAnalysis,'h_rect_prof') 
            
            h_rect_prof = profileAnalysis.h_rect_prof;          
            h_n = str2double(get(hObjs.h_edit_ROI,'String'));
            h_d = str2double(get(hObjs.h_edit_averageWidth,'String'));
                
            pxSzX = imgData.pxSzX;
            
            if round(h_d/pxSzX) > h_n
                
                warndlg(sprintf('max diameter is %d um',floor(h_n*pxSzX)));              
                set(hObjs.h_edit_averageWidth,'String',num2str(floor(h_n*pxSzX)))
                
                return
                                
            end
                       
            if (mod(h_n,2)==0) || (isnan(h_n)) || (~isnumeric(h_n))
               
                warndlg('choose odd number');
                return               
            end   
           
            p_old = getPosition(h_rect_prof);            
            setPosition(h_rect_prof,[p_old(1) p_old(2) p_old(3) h_n])           
            
            plotROIprofile(getPosition(h_rect_prof),mainFig)
            
        else           
            return       
        end
               
end


end





