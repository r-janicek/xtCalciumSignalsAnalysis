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
if isnan(blank)
    errordlg('MEASURE BLANK!')
   return 
end
imgData.blank = blank;
% select type of calculation of normalization 
switch getappdata(mainFig,'analysisType')
    case {'spark recovery ryanodine', 'spark detection'}
        if hObjs.check_bsFit.Value
            typeOfCalc = 2;
        else
            typeOfCalc = 1;
        end
    case 'transients & waves'
        if hObjs.check_doBsFit.Value
            typeOfCalc = 2;
        else
            typeOfCalc = 1;
        end
end
% do pixelwise normalization
% normalize data as: (F-F0)/(F0-blank) ,save F0
% pre-allocation
img_data_n = zeros(size(imgDataXTfluoFN));
img_data_n_raw = zeros(size(imgDataXTfluoFN));
% normalize with fitting or with just mean value division of baseline
% region
switch typeOfCalc
    case 1 % mean value division
        if ~isempty(h_rect_F0)
            % get ROI mask
            roi_maskF0_F0 = createMask(h_rect_F0);
            ROI_pos = h_rect_F0.Position;
            delete(h_rect_F0)
            r_t = find(any(roi_maskF0_F0, 1));
            r_x = find(any(roi_maskF0_F0, 2));
            % get crop of image to calculate baseline values
            crop_F0 = imgDataXTfluoFN(r_x,r_t);
            crop_F0_raw = imgDataXTfluoRN(r_x,r_t);
            % pre-allocation
            F0 = zeros(size(crop_F0,1),1);
            noise_est_px = zeros(size(crop_F0,1),1);
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
                noise_est_px(i,1) = std((F0_line-F0(i,1))./(F0(i,1) - blank), ...
                    [], 'all');
                F0raw(i,1) = mean(F0_line_r);
                % and normalize image to deltaF/F0
                img_data_n(i,:) = (imgDataXTfluoFN(i,:)-F0(i,1)) ./ ...
                    (F0(i,1) - blank);
                img_data_n_raw(i,:) = (imgDataXTfluoRN(i,:)-F0raw(i,1)) ./ ...
                    (F0raw(i,1) - blank);
                clearvars F0_line F0_line_r
            end
            % estimate noise in image
            noise_est = mean(noise_est_px);
        else
            % get detected calcium sparks
            sparkDetection = getappdata(mainFig,'sparkDetection');
            % set areas of detected sparks as NaNs
            mNaN = false(size(imgDataXTfluoFN));
            mNaN(cat(1,sparkDetection.detectedEvents.PixelIdxList)) = true;
            F0 = zeros(size(imgDataXTfluoFN,1),1);
            F0raw = zeros(size(imgDataXTfluoFN,1),1);
            for i = 1:size(imgDataXTfluoFN,1)
                % get baseline estimate as mean of
                % profile without detected events
                F0(i,1) = mean(imgDataXTfluoFN(i,~mNaN(i,:)));
                F0raw(i,1) = mean(imgDataXTfluoRN(i,~mNaN(i,:)));
                % self-ratio
                img_data_n(i,:) = (imgDataXTfluoFN(i,:) - F0(i,1)) ./ ...
                    (F0(i,1) - blank);
                img_data_n_raw(i,:) = (imgDataXTfluoRN(i,:) - F0raw(i,1)) ./ ...
                    (F0raw(i,1) - blank);
            end
            % estimate noise in image
            noise_est = std(img_data_n(~mNaN), [], 'all');
        end

    case 2 % fitting
        % get fit function and baseline mask
        switch getappdata(mainFig,'analysisType')
            case {'spark recovery ryanodine', 'spark detection'}
                typeFit = 'poly1';
                % try to find all events without filtering and
                % applying watershed
                [allEventsMask, ~] = spark_detect_vst(imgDataXTfluoRN, ...
                    10, 3, 8, 1, 5, 0.5);
                mBs = ~imfill(allEventsMask, 'holes');
            case 'transients & waves'
                typeFit = hObjs.popUpMenuBaselineFcn.String{...
                    hObjs.popUpMenuBaselineFcn.Value};
                mBs = imgData.baselineM;
        end
        
        % set up pointer to watch
        mainFig.Pointer = 'watch';
        % waitbar
        hw = waitbar(0,'1','Name','Calculating self-ratio...',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(hw,'canceling',0);
        % pre-allocation
        F0 = zeros(size(imgDataXTfluoFN));
        F0raw = zeros(size(imgDataXTfluoFN));
        for i = 1:size(imgDataXTfluoFN,1)
            % Check for clicked Cancel button
            if getappdata(hw,'canceling')
                break
            end
            % Update waitbar and message
            waitbar(i/size(imgDataXTfluoFN,1), hw, ...
                sprintf('%d %%',round((i/size(imgDataXTfluoFN,1))*100)))
            % calc baseline
            F0(i,:) = fitBaseline(t, imgDataXTfluoFN(i,:), ...
                typeFit, mBs(i,:), 0, []);
            F0raw(i,:) = fitBaseline(t, imgDataXTfluoRN(i,:), ...
                typeFit, mBs(i,:), 0, []);
            % and normalize image to deltaF/F0 
            img_data_n(i,:) = (imgDataXTfluoFN(i,:)-F0(i,:)) ./ ...
                (F0(i,:) - blank);
            img_data_n_raw(i,:) = (imgDataXTfluoRN(i,:)-F0raw(i,:)) ./ ...
                (F0raw(i,:) - blank);
            clearvars F0_line F0_line_r
        end
        delete(hw)
        % set up pointer back
        mainFig.Pointer = 'arrow';
        % estimate noise in image
        noise_est = std(img_data_n(mBs), [], 'all');
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
% crop data in x direction if applicable
if ~isempty(h_rect_F0)
    imgData.imgDataXTfluoR = imgData.imgDataXTfluoR(r_x,:);
    imgData.imgDataXTfluoF = imgData.imgDataXTfluoF(r_x,:);

    imgData.imgDataXT = cellfun(@(x) x(r_x,:), ...
        imgData.imgDataXT, 'UniformOutput',0);
    try
        imgData.imgDataXTtrans = imgData.imgDataXTtrans(r_x,:);
    catch
        imgData.imgDataXTtrans = [];
    end
end
imgData.imgDataXTfluoFN = img_data_n;
imgData.imgDataXTfluoRN = img_data_n_raw;
imgData.F0 = F0;
imgData.F0raw = F0raw;
imgData.z_max_img = z_max_img;
imgData.z_min_img = z_min_img;
imgData.norm_flag = 1;
imgData.stdNoise = noise_est;
% save crop roi position
if ~isempty(h_rect_F0)
    if isfield(imgData,'cropROIpos')
        cropROIpos = imgData.cropROIpos;
        imgData.cropROIpos = [cropROIpos(1)
            cropROIpos(2)+ROI_pos(2)
            cropROIpos(3)
            ROI_pos(4)];
    else
        imgData.cropROIpos = [t(1) ROI_pos(2) t(end)-t(1) ROI_pos(4)];
    end
end
% set up main window
set(hObjs.h_pb_norm, 'String','normalized (F/F0)')
set(hObjs.h_pb_norm, 'Enable','off')
set(hObjs.h_pb_GetBlank, 'Enable','off')
set(hObjs.h_push_BlankROI, 'Enable','off')
set(hObjs.h_edit_Blank, 'Enable','off')
set(hObjs.h_pb_crop, 'Enable','off')
switch getappdata(mainFig,'analysisType')
    case 'spark recovery ryanodine'

    case 'spark detection'
        % % get noise estimate to filter sparks amplitudes
        % % https://www.sciencedirect.com/science/article/pii/S1077314296900600
        % % https://ch.mathworks.com/matlabcentral/fileexchange/36941-fast-noise-estimation-in-images
        % [H W] = size(img_data_n);
        % I = double(img_data_n);
        % % compute sum of absolute values of Laplacian
        % M = [1 -2 1; -2 4 -2; 1 -2 1];
        % Sigma = sum(sum(abs(conv2(I, M))));
        % % scale sigma with proposed coefficients
        % Sigma = Sigma*sqrt(0.5*pi)./(6*(W-2)*(H-2));
        % % noiseEst = std2(img_data_n(~mNaN));
        % imgData.stdNoise = Sigma;
        % set up noise level for spark filtering 2*stdNoise
        hObjs.h_edit_noise.String = num2str(2*imgData.stdNoise);

    case 'transients & waves'
        hObjs.check_doBsFit.Enable = 'off';
        hObjs.popUpMenuBaselineFcn.Enable = 'off';
        hObjs.sld_Bs.Enable = 'off';
end

% save data to main window
setappdata(mainFig, 'imgData', imgData);

end


