function statEvents = findEvents(mainFig,imgRaw)

% find events/sparks in image (not filtered) using method published by:
% Detecting Ca2  sparks on stationary and varying baselines
% Peter Bankhead, C. Norman Scholfield, Tim M. Curtis, and J. Graham McGeown

% img = nonnegative , not filtered image data

% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
analysisType = getappdata(mainFig,'analysisType');
profileAnalysis = getappdata(mainFig,'profileAnalysis');

pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;

%minSpDur = str2double(get(getappdata(mainFig,'h_edit_duration'),'String'));

% get parameters
fFWHM = str2double(get(hObjs.h_edit_fFWHM,'String'));
smoothIter = str2double(get(hObjs.h_edit_smoothIter,'String'));
baseIter = str2double(get(hObjs.h_edit_baseIter,'String'));
tresh = str2double(get(hObjs.h_edit_tresh,'String'));
expFactor = str2double(get(hObjs.h_edit_expFactor,'String'));

watershedTransform = get(hObjs.check_watershed,'Value');
% sld_spDet = hObjs.sld_spDet;
% last pressed pushbutton to detect events
detectEventsPB = getappdata(mainFig,'lastPressedPusbutton'); 
  
% find sparks
[bwImg, ~] = spark_detect_vst(imgRaw, fFWHM, smoothIter, baseIter, 1, tresh, expFactor);

% image of detected events in filtered and normalized image, depends on
% analysis
switch analysisType
    case 'spark detection'
        imgEvents = imgData.imgDataXTfluoFN; % filtered and normalized
        
    case 'spark recovery ryanodine' 
        if strcmp(detectEventsPB.Tag,'findEvents')
            imgEvents = imgData.imgDataXTfluoFN;
        else         
            imgEvents = profileAnalysis.croppedDataWholeAreaRoi;
        end
                   
    otherwise
        imgEvents = imgData.imgDataXTfluoFN;
        
end

% check if there were some events detected
if any(bwImg(:)~=0)
    
    % do watersher and split close merged events if selected
    if watershedTransform
        % spDetS = round(get(hObjs.sld_spDet,'Value')); % in percentage
        % spDetS = get(sld_spDet,'Value');
        % parameters of watershed detection/ also filtering, 
        % wSp in um and dSp in ms
        % filter size is increasing from min values (defined by user) with decreasing sensitivity  
        wSp = max(str2double(hObjs.h_edit_MinWidthSpark.String), 1);
        dSp = max(str2double(hObjs.h_edit_MinDurSpark.String), 10);
        if wSp < pxSzX, wSp = pxSzX; end
        if dSp < pxSzT, dSp = pxSzT; end
        kernelSz_r = 2*ceil(wSp/2/pxSzX)+1;
        kernelSz_c = 2*ceil(dSp/2/pxSzT)+1;
        kernelSz = [kernelSz_r, kernelSz_c];
        % filter with median
        [imgEventsWT, ~] = imgFiltering(imgRaw, pxSzT, pxSzX);
        % calculate self ratio (F/F0) of image pixelwise
        for i = 1:size(imgEventsWT,1)
            % self-ratio
            imgEventsWT(i,:) = imgEventsWT(i,:) ./ ...
                mean(imgEventsWT(i,~bwImg(i,:)));
        end
        % filter with wiener filter
        imgEventsWT = wiener2( ...
            padarray(imgEventsWT, (kernelSz-1)/2, "symmetric"), kernelSz);
        imgEventsWT = imgEventsWT((kernelSz_r-1)/2+1:end-(kernelSz_r-1)/2, ...
            (kernelSz_c-1)/2+1:end-(kernelSz_c-1)/2);
        % filter with gauss
        %imgEventsWT = imgaussfilt(imgEventsWT, 1, "FilterSize",kernelSz);
        % average filter
        imgEventsWT = imboxfilt(imgEventsWT, kernelSz);
        % normalize to [0,1]
        imgEventsWT = (imgEventsWT - min(imgEventsWT(:))) ./ ...
            (max(imgEventsWT(:)) - min(imgEventsWT(:)));
        % imgEventsWT = imgEventsWT ./ max(imgEventsWT(:));
        % % levels from 1 to 50
        % %imgRange = prctile(imgEventsWT(bwImg~=0), [1 99]);
        noiseEst = std(imgEventsWT(~bwImg));
        % nL_min = 1;
        % nL_max = ceil(1/noiseEst);  %3*ceil( (imgRange(2)-imgRange(1)) / noiseEst );
        % convert image to n levels (each level = 2*noise)
        imgEventsWT = imgEventsWT ./ (2*noiseEst);
        % nL = round(((nL_max-nL_min)/(100-1))*spDetS + nL_min - (nL_max-nL_min)/(100-1));
        % % change image to image with nL levels of intensity
        % for r = 1:size(imgEventsWT, 1)
        %     for c = 1:size(imgEventsWT, 2)
        %         if imgEventsWT(r,c) > 0
        %             imgEventsWT(r,c) = round(imgEventsWT(r,c)*nL);
        %         end
        %     end
        % end
        % watershed transform
        L = watershed(-imgEventsWT,8);
        L(L>0)=1;
        bwImg(~logical(L)) = 0;

        % figure
        % tiledlayout(3,1)
        % ax1 = nexttile;
        % imagesc(imgRaw)
        % ax2 = nexttile;
        % imagesc(imgEventsWT)
        % ax3 = nexttile;
        % imagesc(bwImg)
        % linkaxes([ax1, ax2, ax3])
    end
    
    % find all connected components
    CC_event = bwconncomp(bwImg,8);
    
    % calculate properties of all CC (regions)
    statEvents = regionprops(CC_event,imgEvents, 'Centroid','PixelList',...
        'PixelIdxList','SubarrayIdx','PixelValues','MaxIntensity',...
        'WeightedCentroid','BoundingBox','Image','ConvexHull','Area');
    % calculate boundaries of events
    eventsBoundaries = bwboundaries(bwImg,8,'noholes');
    % add them to stat events
    [statEvents.('Boundary')] = eventsBoundaries{:};  
    
    % remove all regions smaller or bigger than specified parameters
    minDur = str2double(hObjs.h_edit_MinDurSpark.String);
    maxDur = str2double(hObjs.h_edit_MaxDurSpark.String);
    minWidth = str2double(hObjs.h_edit_MinWidthSpark.String);
    maxWidth = str2double(hObjs.h_edit_MaxWidthSpark.String);

    m = false(CC_event.NumObjects,1);
    for i=1:CC_event.NumObjects
        
        w = statEvents(i).BoundingBox(3);
        h = statEvents(i).BoundingBox(4);
        
        if w<=ceil(minDur/pxSzT) || h<=ceil(minWidth/pxSzX) || ...
                w>=ceil(maxDur/pxSzT) || h>=ceil(maxWidth/pxSzX)
            m(i,1) = true;
        end
    end
    
    statEvents(m,:) = [];
    
else
    %find all connected components
    CC_event = bwconncomp(bwImg,8);
    statEvents = regionprops(CC_event,imgEvents, 'Centroid','PixelList','PixelIdxList','SubarrayIdx','PixelValues','MaxIntensity',...
        'WeightedCentroid','BoundingBox');
    
end

