function bg_sld_selection(h_O,~)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% handle of figure
h_fig = h_O.Parent;

% get data
if strcmp(h_fig.Name,'sparks analysis')
    
    hObjs = getappdata(h_fig,'hObjs');
    imgData = getappdata(h_fig,'imgData');
    imgSld = hObjs.h_img_sld;
    imgSld_w = str2double(hObjs.h_edit_imgSld_w.String)*1000; % change from s to ms
    h_ax = [hObjs.ax_img,hObjs.ax_img_sparks,hObjs.ax_prof];
    
    hlink1 = hObjs.hlink1;
    hlink2 = hObjs.hlink2;
    
    h_sld_w_e = hObjs.h_edit_imgSld_w;
    h_bg = hObjs.h_bg_sld;
    h_b1 = hObjs.h_b1_sld;
    h_b2 = hObjs.h_b2_sld;
    
else
    hObjs = getappdata(h_fig,'hObjsFit');
    imgData = getappdata(h_fig,'imgData');
    imgSld = hObjs.h_img_sld_fit;
    imgSld_w = str2double(hObjs.h_edit_imgSld_w_fit.String)*1000; % change from s to ms
    h_ax = [hObjs.ax_fit,hObjs.ax_res];
    
    hlink1 = hObjs.hlink;
    try
        hlink2 = hObjs.hlink2;
    catch
    end
    
    h_sld_w_e = hObjs.h_edit_imgSld_w_fit;
    h_bg = hObjs.h_bg_sld_fit;
    h_b1 = hObjs.h_b1_sld_fit;
    h_b2 = hObjs.h_b2_sld_fit;
    
end

if imgSld_w >= imgData.t(end)
    imgSld_w = floor( (imgData.t(end)/2) );
    h_sld_w_e.String = num2str(imgSld_w/1000);
end


switch h_O.Type
    
    case 'uicontrol'
        % show respective part of image
        h_bg.SelectedObject = h_b2;
        set(imgSld,'Value',1,'Max',imgData.t(end)/imgSld_w,'Enable','on')
        % set XLim of axes
        val_l = imgData.t(1);
        val_u = imgData.t(1)+imgSld_w;
        arrayfun(@(x) set(x,'XLim',[val_l val_u]),h_ax)
        
        % remove axes linking, was too slow
        removetarget(hlink1,hlink1.Targets)
        try
            removetarget(hlink2,hlink2.Targets)
        catch
        end
        
        % add listener to slider
        propListener = addlistener(imgSld, 'Value', 'PostSet', @imgShow_slider);
        setappdata(h_fig,'propListener',propListener);
        
        
    otherwise

        switch h_bg.SelectedObject.String
           
            case '#s part'
                % enable slider
                % show respective part of image

                % remove axes linking, was too slow
                removetarget(hlink1,hlink1.Targets)
                try
                    removetarget(hlink2,hlink2.Targets)
                catch
                end
                
                set(imgSld,'Value',1,'Max',imgData.t(end)/imgSld_w,'Enable','on')
                % set XLim of axes
                val_l = imgData.t(1);
                val_u = imgData.t(1)+imgSld_w;
                arrayfun(@(x) set(x,'XLim',[val_l val_u]),h_ax)
                
                % add listener to slider
                propListener = addlistener(imgSld, 'Value', 'PostSet', @imgShow_slider);
                setappdata(h_fig,'propListener',propListener);
 
            otherwise 
                % show full image/profile
            
                % add axes linking
                arrayfun(@(x) addtarget(hlink1,x), h_ax)
                try
                    addtarget(hlink2,hObjs.ax_img_sparks)
                    %addtarget(hlink2,hObjs.ax_img_sparks_2)
                catch
                end
                
                set(imgSld,'Value',1,'Max',ceil(imgData.t(end)/imgSld_w),'Enable','off')
                % set XLim of axes
                val_l = imgData.t(1);
                val_u = imgData.t(end);
                arrayfun(@(x) set(x,'XLim',[val_l val_u]),h_ax)
                
                % remove listener
                delete(getappdata(h_fig,'propListener'))
                
        end
              
end

end

