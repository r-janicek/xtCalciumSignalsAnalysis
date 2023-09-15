function detectEvent(hO,E,analysisFig)

% data
hObjsA = getappdata(analysisFig,'hObjsA');
selectedEvent = getappdata(analysisFig,'selectedEvent');

if ~isempty(hO)
    
    switch hO.Style
        
        case 'pushbutton'
            if isfield(selectedEvent.detectedEvent,'tempMask')
                if ~isempty(selectedEvent.detectedEvent.tempMask)
                    calcMethod = 'usePrevMask';
                else
                    calcMethod = 'calcNewMask';
                end
            else
                calcMethod = 'calcNewMask';
            end
            
        case 'edit'
            calcMethod = 'calcNewMask';
            
    end
    
else
    calcMethod = 'calcNewMask';
end


switch calcMethod
    
    case 'calcNewMask'
        
        %roiName = selectedROIs.roiName{indROI};
        img = selectedEvent.ROIdata.dataROIs.imgFN;
        pxSzT = selectedEvent.pxSzT;
        pxSzX = selectedEvent.pxSzX;
        % t = selectedEvent.ROIdata.dataROIs.t;
        % x = selectedEvent.ROIdata.dataROIs.x;
        
        % get treshold and size filter parameters
        % treshold for image
        F0_prctile = str2double(hObjsA.h_edit_F0_prctile.String);
        trh = str2double(hObjsA.h_edit_treshDet.String);
        fitlSzX = str2double(hObjsA.h_edit_filtSzX.String);  % um
        fitlSzT = str2double(hObjsA.h_edit_filtSzT.String);  % ms
       
        % baseline value
        mI = false(size(img));
        F01 = prctile(img,F0_prctile,2); %mean(img,2);
        % Compute noise estimate from differences along the time dimension
        diffs = diff(img, [], 2);
        % Compute absolute values of required coefficients
        abs_diffs = abs(diffs);
        % Compute noise estimate
        rms = median(abs_diffs,2) ./ sqrt(2) ./ 0.6745;
        % rms = std(img<prctile(img(:),F0_prctile),[],2);
        
        for i=1:size(img,1)
            mL = img(i,:)<(F01(i)+1*rms(i));
            bs = mean( img(i,mL) );
            crit = trh*std( img(i,mL) );
            mI(i,:) = img(i,:)>(bs+crit);
        end
        
        % remove objects with pixels less than: fitlSzX * fitlSzT
        nPix = round( (fitlSzX/pxSzX) * (fitlSzT/pxSzT) );
        BW = bwareaopen(mI,nPix,8);
        % fill holes
        BW = imfill(BW,'holes');
        
        
    case 'usePrevMask'
       img = selectedEvent.ROIdata.dataROIs.imgFN;
       BW = selectedEvent.detectedEvent.tempMask;

end

% find all events and take only biggest
CC = bwconncomp(BW);
stats = regionprops(CC,img,'MaxIntensity','BoundingBox','PixelList','PixelIdxList','SubarrayIdx','Orientation','Area');
fitAreaStat = stats( [stats.Area] == max([stats.Area]) );

% get only biggest event
try
    eventImg = zeros(size(img));
    eventImg(fitAreaStat.PixelIdxList) = img(fitAreaStat.PixelIdxList);
catch 
    eventImg = zeros(size(img));
end

% save
selectedEvent.detectedEvent.mask = BW;
if isfield(selectedEvent.detectedEvent,'tempMask')
    selectedEvent.detectedEvent.tempMask = BW;
end
selectedEvent.detectedEvent.eventImg = eventImg;

setappdata(analysisFig,'selectedEvent',selectedEvent)

% show image of new detected event
try
    h_detImg = findall(hObjsA.ax_detEvent,'Type','image');
    h_detImg.CData = eventImg;
catch 
end
end

