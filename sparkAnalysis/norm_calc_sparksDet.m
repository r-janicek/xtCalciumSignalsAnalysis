function norm_calc_sparksDet(mainFig)
% normalize image, areas of detected sparks are excluded
% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
sparkDetection = getappdata(mainFig,'sparkDetection');

pxSzT = imgData.pxSzT;
t = imgData.t;

imgDataXTfluoFN = imgData.imgDataXTfluoF;
imgDataXTfluoRN = imgData.imgDataXTfluoR;

% get blank 
blank = str2double(get(hObjs.h_edit_Blank,'String'));
if isnan(blank)
    errordlg('MEASURE BLANK!')
   return 
end
imgData.blank = blank;

% set areas of detected sparks as NaNs
mNaN = false(size(imgDataXTfluoFN));
mNaN(cat(1,sparkDetection.detectedEvents.PixelIdxList)) = true;

% do pixelwise normalization
% normalize data as: (F-F0)/(F0-blank) ,save F0
img_data_n = zeros(size(imgDataXTfluoFN));
img_data_n_raw = zeros(size(imgDataXTfluoFN));
F0 = zeros(size(imgDataXTfluoFN,1),1);
F0raw = zeros(size(imgDataXTfluoFN,1),1);

for i = 1:size(imgDataXTfluoFN,1)
    % omnit NaNs
    F0_line = imgDataXTfluoFN(i,~mNaN(i,:));
    F0_line_r = imgDataXTfluoRN(i,~mNaN(i,:));

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
imgData.imgDataXTfluoFN = img_data_n;
imgData.imgDataXTfluoRN = img_data_n_raw;
imgData.F0 = F0;
imgData.F0raw = F0raw;
imgData.z_max_img = z_max_img;
imgData.z_min_img = z_min_img;
imgData.norm_flag = 1;

% get noise estimate
% https://www.sciencedirect.com/science/article/pii/S1077314296900600
% https://ch.mathworks.com/matlabcentral/fileexchange/36941-fast-noise-estimation-in-images
[H W] = size(img_data_n);
I=double(img_data_n);
% compute sum of absolute values of Laplacian
M = [1 -2 1; -2 4 -2; 1 -2 1];
Sigma = sum(sum(abs(conv2(I, M))));
% scale sigma with proposed coefficients
Sigma = Sigma*sqrt(0.5*pi)./(6*(W-2)*(H-2));

% noiseEst = std2(img_data_n(~mNaN));
imgData.stdNoise = Sigma;

% set up noise level for spark filtering 5*stdNoise
hObjs.h_edit_noise.String = num2str(5*Sigma);

% save data
setappdata(mainFig, 'imgData',imgData);

end

