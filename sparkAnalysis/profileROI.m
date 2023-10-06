function profileROI(h_obj,~,mainFig)
% create or change width of ROI to select profile from repetitive sparks
str = get(h_obj,'String');
stl = get(h_obj,'Style');
% get data
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
        if strcmp('whole-cell profile', str)
                            
            if isfield(profileAnalysis,'h_rect_prof')                
                delete(profileAnalysis.h_rect_prof)                           
            end
            keyboard
            h_rect_prof_pos = [0 1 
                size(imgDataXTfluoR,2)*pxSzT 
                size(imgDataXTfluoR,1)];
            plotROIprofile(h_rect_prof_pos,mainFig)
            
        else
            % get width of ROI
            h = str2double(get(hObjs.h_edit_ROI,'String'));
            % get axes limits
            ax_lims = get(ax_img,'YLim');
            
            if isfield(profileAnalysis,'h_rect_prof')                
                delete(profileAnalysis.h_rect_prof)                           
            end
            % ROI
            h_rect_prof = drawrectangle(ax_img, ...
                'Position',[t(1) ax_lims(1) t(end) h], ...
                'Color','r', 'InteractionsAllowed','translate');
            % save ROI
            profileAnalysis.h_rect_prof = h_rect_prof;          
            setappdata(mainFig,'profileAnalysis',profileAnalysis);
            % add interactivity
            addlistener(h_rect_prof, ...
                'ROIMoved',@(src,evnt)plotROIprofile(src,evnt,mainFig));
            set(mainFig, 'WindowKeyPressFcn',...
                {@keyPressProfileROImovement,h_rect_prof,size(imgDataXTfluoR)});
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
           
            p_old = h_rect_prof.Position;     
            h_rect_prof.Position = [p_old(1) p_old(2) p_old(3) h_n];                
            plotROIprofile(h_rect_prof, [], mainFig)        
        else           
            return       
        end         
end
end





