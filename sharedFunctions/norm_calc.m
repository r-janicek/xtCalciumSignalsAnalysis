function norm_calc(st, ~, mainFig, h_rect_F0)
% get data
% style = get(st,'Style');
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
pxSzT = imgData.pxSzT;
t = imgData.t;
imgDataXTfluoFN = imgData.imgDataXTfluoF;
imgDataXTfluoRN = imgData.imgDataXTfluoR;
% get blank
blank = str2double(get(hObjs.h_edit_Blank, 'String'));
imgData.blank = blank;

% get ROI mask
roi_maskF0_F0 = createMask(h_rect_F0);
ROI_pos = h_rect_F0.Position;
delete(h_rect_F0)
r_t = find(any(roi_maskF0_F0, 1));
r_x = find(any(roi_maskF0_F0, 2));
% get crop of image to calculate baseline values
crop_F0 = imgDataXTfluoFN(r_x,r_t);
crop_F0_raw = imgDataXTfluoRN(r_x,r_t);
% select type of calculation of normalization 
switch getappdata(mainFig,'analysisType')
    case 'spark recovery ryanodine'
        typeOfCalc = 1;
    case 'transients & waves'
        if hObjs.check_doBsFit.Value
            typeOfCalc = 2;
        else
            typeOfCalc = 1;
        end
end
% pre-allocation
img_data_n = zeros(size(imgDataXTfluoFN));
img_data_n_raw = zeros(size(imgDataXTfluoFN));
% normalize with fitting or with just mean value division
switch typeOfCalc
    case 1
        % pre-allocation
        F0 = zeros(size(crop_F0,1),1);
        F0raw = zeros(size(crop_F0,1),1);
        % do not take signals higher than k*SD of line for F0 calculations
        % normalize data as: (F-F0)/(F0-blank) ,save F0
        k = 1.5;
        for i = 1:size(crop_F0,1)
            % calc baseline
            F0_line = crop_F0(i,:);
            F0_line = F0_line(F0_line < (k*std(F0_line) + median(F0_line)) );
            F0_line_r = crop_F0_raw(i,:);
            F0_line_r = F0_line_r(F0_line_r < (k*std(F0_line_r) + median(F0_line_r)) );
            F0(i,1) = mean(F0_line);
            F0raw(i,1) = mean(F0_line_r);
            % and normalize image to deltaF/F0
            img_data_n(i,:) = (imgDataXTfluoFN(i,:)-F0(i,1)) ./ ...
                (F0(i,1) - blank);
            img_data_n_raw(i,:) = (imgDataXTfluoRN(i,:)-F0raw(i,1)) ./ ...
                (F0raw(i,1) - blank);
            clearvars F0_line F0_line_r
        end

    case 2 % fitting
        % get fit function
        typeFit = hObjs.popUpMenuBaselineFcn.String{...
            hObjs.popUpMenuBaselineFcn.Value};
        mBs = imgData.baselineM;
        % set up pointer to watch
        mainFig.Pointer = 'watch';
        % waitbar
        hw = waitbar(0,'1','Name','Calculating self-ratio...',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(hw,'canceling',0);
        % pre-allocation
        F0 = zeros(size(imgDataXTfluoFN));
        F0raw = zeros(size(imgDataXTfluoFN));
        for i = 1:size(crop_F0,1)
            % Check for clicked Cancel button
            if getappdata(hw,'canceling')
                break
            end
            % Update waitbar and message
            waitbar(i/size(crop_F0,1), hw, ...
                sprintf('%d %%',round((i/size(crop_F0,1))*100)))
            % calc baseline
            F0(i,:) = fitBaseline(t, imgDataXTfluoFN(i,:), ...
                typeFit, mBs, 0, []);
            F0raw(i,:) = fitBaseline(t, imgDataXTfluoRN(i,:), ...
                typeFit, mBs, 0, []);
            % and normalize image to deltaF/F0 
            img_data_n(i,:) = (imgDataXTfluoFN(i,:)-F0(i,:)) ./ ...
                (F0(i,:) - blank);
            img_data_n_raw(i,:) = (imgDataXTfluoRN(i,:)-F0raw(i,1)) ./ ...
                (F0raw(i,1) - blank);
            clearvars F0_line F0_line_r
        end
        delete(hw)
        % set up pointer back
        mainFig.Pointer = 'arrow';    
end

% set NaN as zero
img_data_n(isnan(img_data_n)) = 0;
img_data_n_raw(isnan(img_data_n_raw)) = 0;

mm = sort(max(img_data_n,[],2));
z_max_img = mm(round(0.90*length(max(img_data_n,[],2))));
z_min_img = min(min(img_data_n));
%z_max_img = prctile(img_data_n,99,'all');
%z_min_img = prctile(img_data_n,1,'all');

% plot normalized image 
ax_img = hObjs.ax_img;
whole_img_h = findobj(ax_img, 'Type','Image');
set(whole_img_h, 'CData',img_data_n, ...
    'XData',[0 max(t)], ...
    'YData',[1 size(img_data_n,1)])
%set(ax_img,'CLim',prctile(img_data_n(:),[0.1 99.9]))
set(ax_img, 'CLim',[z_min_img z_max_img])
set(ax_img, 'XLim',[0 max(t)], 'YLim',[1 size(img_data_n,1)])
set(ax_img, 'XTick',[])
set(get(ax_img, 'Ylabel'), 'String','x (pixels) [filtered & normalized]')

% plot normalized time profile
ax_prof = hObjs.ax_prof;
plot(t, mean(img_data_n,1), 'Parent',ax_prof)
set(ax_prof, 'Xlim',[t(1) t(end)], ...
    'YLim', getAxisLimits(mean(img_data_n, 1), 5),...
    'FontSize',14)
set(get(ax_prof,'Xlabel'), 'String','t (ms)', 'FontWeight','bold')
set(get(ax_prof,'Ylabel'), 'String','fluorescence (deltaF/F0)', ...
    'FontWeight','bold')

% save data 
% crop data in x direction
imgData.imgDataXTfluoR = imgData.imgDataXTfluoR(r_x,:);
imgData.imgDataXTfluoF = imgData.imgDataXTfluoF(r_x,:);
try
    imgData.imgDataXTtrans = imgData.imgDataXTtrans(r_x,:);
catch
    imgData.imgDataXTtrans = [];
end
imgData.imgDataXT = cellfun(@(x) x(r_x,:), ...
    imgData.imgDataXT, 'UniformOutput',0);

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
    imgData.cropROIpos = [cropROIpos(1) 
                          cropROIpos(2)+ROI_pos(2) 
                          cropROIpos(3) 
                          ROI_pos(4)];
else
    imgData.cropROIpos = [t(1) ROI_pos(2) t(end)-t(1) ROI_pos(4)];
end

setappdata(mainFig, 'imgData', imgData);

set(hObjs.h_pb_norm, 'String','normalized (F/F0)')
set(hObjs.h_pb_norm, 'Enable','off')

set(hObjs.h_pb_GetBlank, 'Enable','off')
set(hObjs.h_push_BlankROI, 'Enable','off')
set(hObjs.h_edit_Blank, 'Enable','off')

set(hObjs.h_edit_Blank, 'Enable','off')
switch getappdata(mainFig,'analysisType')
    case 'spark recovery ryanodine'

    case 'transients & waves'
        hObjs.check_doBsFit.Enable = 'off';
        hObjs.popUpMenuBaselineFcn.Enable = 'off';
        hObjs.sld_Bs.Enable = 'off';
end

end


