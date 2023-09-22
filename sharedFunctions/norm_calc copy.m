function norm_calc(st,~,mainFig,h_rect_F0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% get data
% style = get(st,'Style');
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

pxSzT = imgData.pxSzT;
t = imgData.t;

imgDataXTfluoFN = imgData.imgDataXTfluoF;
imgDataXTfluoRN = imgData.imgDataXTfluoR;

% get blank 
blank = str2double(get(hObjs.h_edit_Blank,'String'));
imgData.blank = blank;

% get ROI position
ROI_pos = getPosition(h_rect_F0);
pos = [ROI_pos(1)/pxSzT ROI_pos(2) ROI_pos(3)/pxSzT+0.5 ROI_pos(4)];
pos(pos<0)=0;
delete(h_rect_F0)

[r_t(:,1),r_x(:,1)] = rect2ind(pos);

if length(r_t)>size(imgDataXTfluoFN,2)
    r_t = r_t(1:size(imgDataXTfluoFN,2),1);
end

if length(r_x)>size(imgDataXTfluoFN,1)
    r_x = r_x(1:size(imgDataXTfluoFN,1),1);
end

r_t(r_t > size(imgDataXTfluoFN,2))=[];
r_x(r_x > size(imgDataXTfluoFN,1))=[];

crop_F0 = imgDataXTfluoFN(r_x,r_t);
imgDataXTfluoFN = imgDataXTfluoFN(r_x,:);

crop_F0_raw = imgDataXTfluoRN(r_x,r_t);
imgDataXTfluoRN = imgDataXTfluoRN(r_x,:);

img_data_n = zeros(size(imgDataXTfluoFN));
img_data_n_raw = zeros(size(imgDataXTfluoFN));
F0 = zeros(size(crop_F0,1),1);
F0raw = zeros(size(crop_F0,1),1);
% do not take signals higher than k*SD of line for F0 calculations
% normalize data as: (F-F0)/(F0-blank) ,save F0
k = 1.5;
for i = 1:size(crop_F0,1)
  
    F0_line = crop_F0(i,:);
    F0_line = F0_line(F0_line < (k*std(F0_line) + median(F0_line)) );
    
    F0_line_r = crop_F0_raw(i,:);
    F0_line_r = F0_line_r(F0_line_r < (k*std(F0_line_r) + median(F0_line_r)) );
  
    F0(i,1) = mean(F0_line);
    F0raw(i,1) = mean(F0_line_r);
        
    img_data_n(i,:) = (imgDataXTfluoFN(i,:)-F0(i,1)) ./ (F0(i,1) - blank);
    img_data_n_raw(i,:) = (imgDataXTfluoRN(i,:)-F0raw(i,1)) ./ (F0raw(i,1) - blank);
        
    clearvars F0_line F0_line_r 
end

% set NaN as zero
img_data_n(isnan(img_data_n))=0;
img_data_n_raw(isnan(img_data_n_raw))=0;

mm = sort(max(img_data_n,[],2));
z_max_img = mm(round(0.90*length(max(img_data_n,[],2))));
z_min_img = min(min(img_data_n));

% plot normalized image 
ax_img = hObjs.ax_img;
whole_img_h = findobj(ax_img,'Type','Image');
set(whole_img_h,'CData',img_data_n,'XData',[0 max(t)],'YData',[1 size(img_data_n,1)])
%set(ax_img,'CLim',prctile(img_data_n(:),[0.1 99.9]))
set(ax_img,'CLim',[z_min_img z_max_img])
set(ax_img,'XLim',[0 max(t)],'YLim',[1 size(img_data_n,1)])

set(ax_img,'XTick',[])
set(get(ax_img,'Ylabel'),'String','x (pixels) [filtered & normalized]')

% plot normalized time profile
ax_prof = hObjs.ax_prof;
plot(t,mean(img_data_n,1),'Parent',ax_prof)
set(ax_prof,'Xlim',[t(1) t(end)],'YLim',[min(mean(img_data_n,1)) max(mean(img_data_n,1))],...
    'FontSize',14)
set(get(ax_prof,'Xlabel'),'String','t (ms)','FontWeight','bold')
set(get(ax_prof,'Ylabel'),'String','Fluorescence (F/F0)','FontWeight','bold')

% save data 
% crop data in x direction
imgData.imgDataXTfluoR = imgData.imgDataXTfluoR(r_x,:);
imgData.imgDataXTfluoF = imgData.imgDataXTfluoF(r_x,:);

try
    imgData.imgDataXTtrans = imgData.imgDataXTtrans(r_x,:);
catch
    imgData.imgDataXTtrans = [];
end

try
    imgData.imgDataXT = cellfun(@(x) x(r_x,:),imgData.imgDataXT,'UniformOutput',0);
catch 
    imgData.imgDataXT = imgData.imgDataXT(r_x,:);
end

imgData.imgDataXTfluoFN = img_data_n;
imgData.imgDataXTfluoRN = img_data_n_raw;
imgData.F0 = F0;
imgData.F0raw = F0raw;
imgData.z_max_img = z_max_img;
imgData.z_min_img = z_min_img;
imgData.norm_flag = 1;

% save crop roi position
if isfield(imgData,'cropROIpos')
    cropROIpos = imgData.cropROIpos;   
    imgData.cropROIpos = [cropROIpos(1) cropROIpos(2)+ROI_pos(2) cropROIpos(3) ROI_pos(4)];
else
    imgData.cropROIpos = [t(1) ROI_pos(2) t(end)-t(1) ROI_pos(4)];
end

% correct positions of lines or rectangles in cropped image showing photolytic pulses
lines_TPP = findobj(ax_img,'Type','Line','Tag','photolysis');
rect_TPP = findobj(ax_img,'Type','Rectangle','Tag','photolysis');

if ~isempty(lines_TPP)
    arrayfun(@(x) set(x,'YData',[1 size(img_data_n,1)]), lines_TPP, 'UniformOutput', false)
end

if ~isempty(rect_TPP)   
    arrayfun(@(x) set(x,'Position',x.Position - [0  ROI_pos(2) 0 0]), rect_TPP, 'UniformOutput', false)   
end

%stdNoise = std2(img_data_n((img_data_n<(1+2*std2(img_data_n)))&(img_data_n>(1-2*std2(img_data_n)))));

setappdata(mainFig,'imgData',imgData);

set(hObjs.h_pb_norm,'String','normalized (F/F0)')
set(hObjs.h_pb_norm,'Enable','off')

set(hObjs.h_pb_GetBlank,'Enable','off')
set(hObjs.h_push_BlankROI,'Enable','off')
set(hObjs.h_edit_Blank,'Enable','off')

set(hObjs.h_edit_Blank,'Enable','off')

end


