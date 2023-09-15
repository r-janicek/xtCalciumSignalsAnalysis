function mouseProfileRoiMovement(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mainFig = varargin{1,1};
currObj = get(mainFig,'CurrentObject');
hObjs = getappdata(mainFig,'hObjs');
profileAnalysis = getappdata(mainFig,'profileAnalysis');

if ~isempty(currObj)   
    if strcmp(currObj.Type,'patch')

        ax_img = hObjs.ax_img;
        h_rect_prof = profileAnalysis.h_rect_prof;       
        h_rect_prof_pos = getPosition(h_rect_prof);
              
        plotROIprofile(h_rect_prof_pos,mainFig)
       
    end
end
