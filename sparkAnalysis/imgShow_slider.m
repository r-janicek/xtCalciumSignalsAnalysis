function imgShow_slider(varargin)  %imgSld,~
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
try
% handle of slider
imgSld = varargin{1,2}.AffectedObject;
catch 
   keyboard 
end
% handle of figure
h_fig = imgSld.Parent;

% get data
if strcmp(h_fig.Name,'sparks analysis')
    hObjs = getappdata(h_fig,'hObjs');
    imgSld_w = str2double(hObjs.h_edit_imgSld_w.String)*1000; % change from s to ms
    h_ax = [hObjs.ax_img,hObjs.ax_img_sparks,hObjs.ax_prof];
    
else 
    hObjs = getappdata(h_fig,'hObjsFit');
    imgSld_w = str2double(hObjs.h_edit_imgSld_w_fit.String)*1000; % change from s to ms
    h_ax = [hObjs.ax_fit,hObjs.ax_res];
    
end

% update axes
if strcmp(imgSld.Enable,'on')
    % set x limits of axes
    val_l = imgSld.Value*imgSld_w-imgSld_w;
    val_u = imgSld.Value*imgSld_w;
    
    arrayfun(@(x) set(x,'XLim',[val_l val_u]),h_ax)
   
end

end



