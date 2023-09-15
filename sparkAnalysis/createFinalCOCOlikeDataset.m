function createFinalCOCOlikeDataset(~, ~, mainFig)
% create final datasets from already processed images
% also create json file for each dataset with annotations
% create 3 datasets: 60% train, 20% val, 20% test

%% open folder from which to take data of images (.mat files) to create dataset
% load whole selected directory
hObjs = getappdata(mainFig,'hObjs');

folderPath = uigetdir('/Users/rjanicek/Desktop/','Analyze all images from the folder');
if folderPath == 0
    return
end

% change pointer to watch
set(mainFig,'Pointer','watch')
drawnow


%% create folders for data and delete old ones
% create new folder with cutted/splitted images
folderPath_splittedImgs = fullfile(folderPath,'splittedImgs');
if exist(folderPath_splittedImgs, 'dir')
    rmdir(folderPath_splittedImgs, 's')
end
mkdir(folderPath_splittedImgs)

% create folders for train, val and test
folderPath_train = fullfile(folderPath,'train');
folderPath_val = fullfile(folderPath,'val');
folderPath_test = fullfile(folderPath,'test');

if exist(folderPath_train, 'dir')
    rmdir(folderPath_train, 's')
end
if exist(folderPath_val, 'dir')
    rmdir(folderPath_val, 's')
end
if exist(folderPath_test, 'dir')
    rmdir(folderPath_test, 's')
end

mkdir(folderPath_train)
mkdir(folderPath_val)
mkdir(folderPath_test)

% get list of paths of all images data 
fileList = getAllFilesFromFolder(folderPath);
processedImgsPaths = fileList(~cellfun(@isempty,regexp(fileList,'.mat')));


%% split images and masks into parts with maximum defined size
% cut images so max dimension lenght is 256 px for imgs with sparks and 1024 for imgs with waves and transients  
finalImgSize = str2double(hObjs.h_maxImgDimSz_edit.String);
removeCategories = [2,3,6];
maxImgSz_wave_transient = 1024;

% add waitbar 
hw = waitbar(0,'processing image #1','Name','splitting images into parts ...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(hw,'canceling',0)

splittedImgsPaths = {};
getAllCategoriesInImgs = [];
imgNum = 1;
for i = 1:numel(processedImgsPaths)
    
    % load image data 
    load(processedImgsPaths{i});
    img  = savedImgData.img;
    mask = savedImgData.maskOfEvents;
    pxSzT = savedImgData.pxSzT;
    pxSzX = savedImgData.pxSzX;
    
    % check if transients-5 are present, 
    % if yes cut them out and save them, so they are not cutted in the
    % middle, transient propagates in whole cell, so cut only in t
    % dimension
    if any(ismember(mask(:), 5))
        masks_all_transients_in_t = any(mask==5,1);
        transients_BB = regionprops(mask==5, 'BoundingBox');
        % take the middle of transient and take the rectangle around it in
        % that way, so there is no cut in any other transient mask
        for t = 1:numel(transients_BB)
            tr_BB = transients_BB(t).BoundingBox;
            % middle point of transient
            m_tr = tr_BB(1)+tr_BB(3)/2;
            % mask of suggested cut out for transient
            try
                mask_transient_in_t = false(size(masks_all_transients_in_t));
                masks_allOthers_transients_in_t  = masks_all_transients_in_t;
                masks_allOthers_transients_in_t(round(tr_BB(1)):round(tr_BB(1)+tr_BB(3))) = false;
            catch
                keyboard
            end
            % cut out transient image
            tr_s = round(m_tr-maxImgSz_wave_transient/2);
            if tr_s<1, tr_s = 1; end
            tr_e = round(m_tr+maxImgSz_wave_transient/2);
            if tr_e>size(mask,2), tr_e = size(mask,2); end
            mask_transient_in_t(tr_s:tr_e) = true;  
            % check for overlaps with other transients
            mask_transient_in_t = mask_transient_in_t & ...
                (~masks_allOthers_transients_in_t);
            % adjusted start and end
            tr_s = find(mask_transient_in_t,1,'first');
            tr_e = find(mask_transient_in_t,1,'last');
            
            % save cutted images of transients 
            % new name for part of image
            newImgName = sprintf('imgData_%d.mat',imgNum);
            
            % set up path of created folder
            pathNewImg = fullfile(folderPath_splittedImgs, newImgName);
            
            % resize transient images to selected size
            [tr_img, tr_mask] = resizeImgToSpecificSize( ...
                img(:, tr_s:tr_e), mask(:, tr_s:tr_e), finalImgSize);

            % check if by any chance transient is in begining
            mask_sum = sum(tr_mask, 1);
            if mask_sum(1)>1
                continue
            end
            
            % save image and mask as structure
            splittedImgData = struct( ...
                'img', tr_img, ...
                'maskOfEvents', tr_mask, ...
                'pxSzX', pxSzX, ...
                'pxSzT', pxSzT );
            save(pathNewImg, 'splittedImgData')
            
            % save paths of splitted images
            splittedImgsPaths = [splittedImgsPaths; pathNewImg];
            
            % get all categories in dataset, here only transients-5
            getAllCategoriesInImgs = [getAllCategoriesInImgs; 5];
            
            imgNum = imgNum+1;
        end
        
        % update img and mask, so to process only parts where are no transients
        img = img(:,tr_e:end);
        mask = mask(:,tr_e:end);
    end
    
    % once saved, remove transients and other unwanted categories from mask
    mask(ismember(mask,[removeCategories, 5])) = 0;
       
    switch join(string(unique(mask)),'')
        case '0'
            % no events
            continue
            
        case '01'
            % sparks only
            maxImgSz = 256;
            [imgNum, splittedImgsPaths, getAllCategoriesInImgs] = ...
                splitImageAndMaskIntoParts(img, mask, maxImgSz, imgNum, ...
                splittedImgsPaths, getAllCategoriesInImgs, ...
                folderPath_splittedImgs, pxSzX, pxSzT, finalImgSize, ...
                removeCategories);

        case '04'
            % waves only
            maxImgSz = 1024;
            [imgNum, splittedImgsPaths, getAllCategoriesInImgs] = ...
                splitImageAndMaskIntoParts(img, mask, maxImgSz, imgNum, ...
                splittedImgsPaths, getAllCategoriesInImgs, ...
                folderPath_splittedImgs, pxSzX, pxSzT, finalImgSize, ...
                removeCategories);
            
        case '014'
            % waves and sparks
            maxImgSz = 256;
            [imgNum, splittedImgsPaths, getAllCategoriesInImgs] = ...
                splitImageAndMaskIntoParts(img, mask, maxImgSz, imgNum, ...
                splittedImgsPaths, getAllCategoriesInImgs, ...
                folderPath_splittedImgs, pxSzX, pxSzT, finalImgSize, ...
                removeCategories);
            maxImgSz = 1024;
            [imgNum, splittedImgsPaths, getAllCategoriesInImgs] = ...
                splitImageAndMaskIntoParts(img, mask, maxImgSz, imgNum, ...
                splittedImgsPaths, getAllCategoriesInImgs, ...
                folderPath_splittedImgs, pxSzX, pxSzT, finalImgSize, ...
                removeCategories);
                
    end

    % update waitbar's message field
    waitbar( i/numel(processedImgsPaths),hw,...
        sprintf( 'processing image #%d', i) )
    
    % Check for clicked Cancel button
    if getappdata(hw,'canceling')
        % delete waitbar
        delete(hw)
        break
    end
    
end

% delete waitbar
delete(hw), clearvars hw

% all present categories of events
allCategoriesInImgs = unique(getAllCategoriesInImgs);
allCategoriesInImgs(allCategoriesInImgs==0) = [];


%% create 3 datasets: 60% train, 20% val, 10% test + respective json files
% for annotations

% create structure for json file. COCO-like
% info dataset
infoDataset_json = struct( ...
    'description', 'calciumSignals', ...
    'url', '', ...
    'version', 1.0, ...
    'year', '2021', ...
    'contributor', 'Rado Janicek', ...
    'date_created', datestr(today('datetime')) );

% licenses for dataset
licenseDataset_json = table( ...
    {''}, 1, {'license'} ,...
    'VariableNames',{'url', 'id', 'name'});

% categories
categoriesInDataset_json = table( ...
        {'calciumSignals';'calciumSignals';'calciumSignals'; ...
        'calciumSignals';'calciumSignals';'calciumSignals'},...
        [1;2;3;4;5;6], ...
        {'spark'; 'longLastingSpark'; 'miniWave'; ...
        'wave'; 'transient'; 'caffeineTransient'} ,...
        'VariableNames',{'supercategory', 'id', 'name'} );

categoriesInDataset_json( ...
    ~ismember(categoriesInDataset_json.id, allCategoriesInImgs),:) = [];

% change categories id to 1,2,3,...
oldCategories = table( ... 
    categoriesInDataset_json.id, ...
    (1:height(categoriesInDataset_json))', ...
    'VariableNames', {'oldCatID', 'newCatID'});
categoriesInDataset_json.id = (1:height(categoriesInDataset_json))';

% initialize structures
trainImgs = struct( ...
    'id', {}, ...
    'license', {}, ...
    'coco_url', {}, ...
    'flickr_url', {}, ...
    'width', {}, ...
    'height', {}, ...
    'file_name', {}, ...
    'date_captured', {} );

valImgs = struct( ...
    'id', {}, ...
    'license', {}, ...
    'coco_url', {}, ...
    'flickr_url', {}, ...
    'width', {}, ...
    'height', {}, ...
    'file_name', {}, ...
    'date_captured', {} );

testImgs = struct( ...
    'id', {}, ...
    'license', {}, ...
    'coco_url', {}, ...
    'flickr_url', {}, ...
    'width', {}, ...
    'height', {}, ...
    'file_name', {}, ...
    'date_captured', {} );

annotationsInTrainImgs = struct( ...
    'id', {}, ...
    'category_id', {}, ...
    'iscrowd', {}, ...
    'segmentation', {}, ... % [x1,y1,x2,y2,...]
    'image_id', {}, ...
    'area', {}, ... % number of pixels
    'bbox', {} );  %  boundig box [x,y,w,h] [0 0] is left upper corner

annotationsInValImgs = struct( ...
    'id', {}, ...
    'category_id', {}, ...
    'iscrowd', {}, ...
    'segmentation', {}, ... % [x1,y1,x2,y2,...]
    'image_id', {}, ...
    'area', {}, ... % number of pixels
    'bbox', {} );  %  boundig box [x,y,w,h] [0 0] is left upper corner

annotationsInTestImgs = struct( ...
    'id', {}, ...
    'category_id', {}, ...
    'iscrowd', {}, ...
    'segmentation', {}, ... % [x1,y1,x2,y2,...]
    'image_id', {}, ...
    'area', {}, ... % number of pixels
    'bbox', {} );  %  boundig box [x,y,w,h] [0 0] is left upper corner

% annotation id, starts at 10000
annotID = 10000;

%% main loop, get description from saved splitted images and analyze images for sparks
% add waitbar 
hw = waitbar(0,'0 %; 0 s','Name','creating train, test and val data ...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(hw,'canceling',0)

allSparks = table();
tIter = nan(numel(splittedImgsPaths),1);
% shuffle imgs
idx_imgsPaths = randperm(numel(splittedImgsPaths));
for i = 1:numel(idx_imgsPaths)
    
    % get time per iteration
    tStartIter = tic;
    
    idx = idx_imgsPaths(i);
    
    % load data of splitted imgages
    load(splittedImgsPaths{idx})
    
    imgSplitted = splittedImgData.img;
    maskSplitted = splittedImgData.maskOfEvents;
    
    if numel(unique(maskSplitted))<2
        keyboard
    end
    
    
    % get properties of regions in image
    %find all connected components
    CC_event = bwconncomp(maskSplitted, 8);
    
    %calculate properties of all CC (regions)
    statEvents = regionprops(CC_event, maskSplitted, 'PixelList',...
        'PixelIdxList', 'SubarrayIdx', 'MaxIntensity',...
        'BoundingBox', 'ConvexHull', 'Area');
    % use boundary from find events, not convex hull property
    % calculate boundaries of events
    eventsBoundaries = bwboundaries(maskSplitted,8,'noholes');
    % add them to stat events
    [statEvents.('Boundary')] = eventsBoundaries{:};  
    

    if ~any(ismember(oldCategories.oldCatID, [statEvents.MaxIntensity] ))
        keyboard
    end
    
    
%     newCatID = oldCategories.newCatID( ...
%             oldCategories.oldCatID == statEvents.MaxIntensity );
    
    % skip empty images
    if isempty(statEvents)
        keyboard
        continue
    end
    
    % image structure
    img_json = struct( ...
        'id', idx, ...
        'license', licenseDataset_json.id, ...
        'coco_url', '', ...
        'flickr_url', '', ...
        'width', size(imgSplitted,2), ...
        'height', size(imgSplitted,1), ...
        'file_name', sprintf('img_%d.png',idx), ...
        'date_captured', '');
    
    % create structure for all anotations in image
    allAnnotationsInImg_json = struct( ...
        'id', {}, ...
        'category_id', {}, ...
        'iscrowd', {}, ...
        'segmentation', {}, ... % [x1,y1,x2,y2,...]
        'image_id', {}, ...
        'area', {}, ... % number of pixels
        'bbox', {} );  %  boundig box [x,y,w,h] [0 0] is left upper corner
    for j = 1:numel(statEvents)
        
        %eventBorder = statEvents(j).ConvexHull;
        eventBorder = [ statEvents(j).Boundary(:,2), ...
                        statEvents(j).Boundary(:,1) ];
        bbox = statEvents(j).BoundingBox;
        % in python idx starts at 0
        bbox = [ bbox(1)-1, bbox(2)-1, bbox(3), bbox(4) ];
        eventBorder = eventBorder-1;
        eventBorder = eventBorder';
        eventBorder = eventBorder(:);
        eventBorder = eventBorder';
        
        newCatID = oldCategories.newCatID( ...
            oldCategories.oldCatID == statEvents(j).MaxIntensity );
        
        if isempty(newCatID)
            continue
        end
        
        singleAnnotationsInImg = struct( ...
            'id', annotID, ...
            'category_id', newCatID, ... % statEvents(j).MaxIntensity, ...
            'iscrowd', 0, ...
            'segmentation', {{eventBorder}}, ... % [x1,y1,x2,y2,...]
            'image_id', idx, ...
            'area', statEvents(j).Area, ... % number of pixels
            'bbox', bbox );  %  boundig box [x,y,w,h] [0 0] is left upper corner
        
        allAnnotationsInImg_json = ...
            [allAnnotationsInImg_json; singleAnnotationsInImg];
        % add 1 to annotaion ID
        annotID = annotID+1;
    end

    % first 70% of images as train, 20% as validation and 10% as test dataset
    if i <= round(numel(idx_imgsPaths)*0.7)
        % save image as png
        imwrite(imgSplitted, ...
            fullfile(folderPath_train,img_json.file_name), ...
            'BitDepth',16);
        trainImgs = [trainImgs; img_json];
        annotationsInTrainImgs = ...
            [annotationsInTrainImgs; allAnnotationsInImg_json];
        
    elseif ( i > round(numel(idx_imgsPaths)*0.7) && ...
            i <= round(numel(idx_imgsPaths)*0.9) )
        % save image as png
        imwrite(imgSplitted, ...
            fullfile(folderPath_val,img_json.file_name), ...
            'BitDepth',16);
        valImgs = [valImgs ; img_json];
        annotationsInValImgs = ...
            [annotationsInValImgs; allAnnotationsInImg_json];
        
    else
        % save image as png
        imwrite(imgSplitted, ...
            fullfile(folderPath_test,img_json.file_name), ...
            'BitDepth',16);
        testImgs = [testImgs; img_json];
        annotationsInTestImgs = ...
            [annotationsInTestImgs; allAnnotationsInImg_json];
        
    end
    
    % get time per iteration, in seconds
    tIter(i) = toc(tStartIter);
    estimTimeOfFinish = mean(tIter,'omitnan')*(numel(idx_imgsPaths)-i);
    % in hours
    estimTimeOfFinish_hour = [floor(estimTimeOfFinish/3600), ...
        round(mod(estimTimeOfFinish,3600)/60)];
    % in minutes
    estimTimeOfFinish_min = floor(estimTimeOfFinish/60);
    % in seconds
    estimTimeOfFinish_sec = round(estimTimeOfFinish); % round(mod(estimTimeOfFinish,60));
    
    if estimTimeOfFinish_hour(1) >= 1
        % Report current estimate in the waitbar's message field
        waitbar( i/numel(idx_imgsPaths),hw,...
            sprintf( '%d %%; remaining time: %d hours and %d minutes', ...
            round(100*i/numel(idx_imgsPaths)),...
            estimTimeOfFinish_hour(1),estimTimeOfFinish_hour(2) ) )
    else
        if estimTimeOfFinish_min > 1
            % Report current estimate in the waitbar's message field
            waitbar( i/numel(idx_imgsPaths),hw,...
                sprintf( '%d %%; remaining time: %d minutes', ...
                round(100*i/numel(idx_imgsPaths)),...
                estimTimeOfFinish_min ) )
        else
            % Report current estimate in the waitbar's message field
            waitbar( i/numel(idx_imgsPaths),hw,...
                sprintf( '%d %%; remaining time: %d seconds', ...
                round(100*i/numel(idx_imgsPaths)),...
                estimTimeOfFinish_sec ) )
        end
    end
    
    % Check for clicked Cancel button
    if getappdata(hw,'canceling')
        % delete waitbar
        delete(hw)
        break
    end
    
end

% delete waitbar
delete(hw), clearvars hw


%% create and save final json file with annotations
% with name of folder and save it into folder
final_json_file_train = jsonencode( struct( ...
    'info', infoDataset_json, ...
    'licenses', licenseDataset_json, ...
    'images', trainImgs, ...
    'annotations', annotationsInTrainImgs, ...
    'categories', categoriesInDataset_json ) );

final_json_file_val = jsonencode( struct( ...
    'info', infoDataset_json, ...
    'licenses', licenseDataset_json, ...
    'images', valImgs, ...
    'annotations', annotationsInValImgs, ...
    'categories', categoriesInDataset_json ) );

final_json_file_test = jsonencode( struct( ...
    'info', infoDataset_json, ...
    'licenses', licenseDataset_json, ...
    'images', testImgs, ...
    'annotations', annotationsInTestImgs, ...
    'categories', categoriesInDataset_json ) );

% save json files
if exist(fullfile(folderPath_train,'annotations.json'), 'file')
    delete(fullfile(folderPath_train,'annotations.json'))
end
fid_train = fopen( ...
    fullfile(folderPath_train,'annotations.json'), 'w');
fwrite(fid_train, final_json_file_train, 'char');
fclose(fid_train);

if exist(fullfile(folderPath_val,'annotations.json'), 'file')
    delete(fullfile(folderPath_val,'annotations.json'))
end
fid_val = fopen( ...
    fullfile(folderPath_val,'annotations.json'), 'w');
fwrite(fid_val, final_json_file_val, 'char');
fclose(fid_val);

if exist(fullfile(folderPath_test,'annotations.json'), 'file')
    delete(fullfile(folderPath_test,'annotations.json'))
end
fid_test = fopen( ...
    fullfile(folderPath_test,'annotations.json'), 'w');
fwrite(fid_test, final_json_file_test, 'char');
fclose(fid_test);

% remove folder with splitted images data
rmdir(folderPath_splittedImgs, 's')

% put back pointer to arrow
set(mainFig,'Pointer','arrow')
drawnow  

end


%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
function listOfFiles = getAllFilesFromFolder(folderPath)

    % get content of folder
    listing = dir(folderPath);

    % list of files
    listOfFiles = {listing(~[listing.isdir]).name}';
    % add path to file
    listOfFiles = cellfun(@(x) fullfile(folderPath,x),listOfFiles,'UniformOutput',0);
    
    % look for subdirectories
    subDirs = {listing([listing.isdir]).name};
    % do not take . or ..
    validSubDirs = subDirs(~ismember(subDirs,{'.','..'}));
    
    % loop through folder to look for subfolders
    for subDir = validSubDirs
        listOfFiles = [listOfFiles;getAllFilesFromFolder(fullfile(folderPath,subDir{1}))];
    end
 
end


function [imgNum, splittedImgsPaths, getAllCategoriesInImgs] = ...
        splitImageAndMaskIntoParts(img, mask, maxImgSz, imgNum, ...
        splittedImgsPaths, getAllCategoriesInImgs, ...
        folderPath_splittedImgs, pxSzX, pxSzT, finalImgSize, ...
        removeCategories)
    
    % check size of image, if necessary to split it
    if any( size(img) > maxImgSz )
        % split image into parts
        a = size(img, 1);
        b = size(img, 2);
        
        nParts_a = floor(a/maxImgSz);
        rem_a = rem(a, maxImgSz);
        if nParts_a == 0
            nParts_a = 1;
            a_sz = ones(1,nParts_a)*rem_a;
        else
            a_sz = [ones(1,nParts_a)*maxImgSz, rem_a];
        end
        
        nParts_b = floor(b/maxImgSz);
        rem_b = rem(b, maxImgSz);
        if nParts_b == 0
            nParts_b = 1;
            b_sz = ones(1,nParts_b)*rem_b;
        else
            b_sz = [ones(1,nParts_b)*maxImgSz, rem_b];
        end
        
        imgParts = mat2cell(img, a_sz, b_sz);
        maskParts = mat2cell(mask, a_sz, b_sz);
        
    else
        imgParts = {img};
        maskParts = {mask};
    end
    
    % do not save very small parts
    m_sz_imgParts = cellfun(@(x) min(size(x)),imgParts) < 100;
    imgParts = imgParts(~m_sz_imgParts);
    maskParts = maskParts(~m_sz_imgParts);
    
%     % do not take specific categories if defined
%     m_removeCategories = cellfun( ...
%         @(x) any( ismember( unique(x(:)), removeCategories ) ), maskParts,...
%         'UniformOutput', true);
    
%     imgParts = imgParts(~m_removeCategories);
%     maskParts = maskParts(~m_removeCategories);
    
    % get all categories in dataset
    ctgInImg = [];
    for ctg =1:numel(maskParts)
        ctgInImg = [ctgInImg; unique(maskParts{ctg})];
    end
    getAllCategoriesInImgs = [getAllCategoriesInImgs; ctgInImg];
    
    for p = 1:numel(imgParts)
        
        % resize to 256 pixels
        [imgParts{p}, maskParts{p}] = resizeImgToSpecificSize( ...
            imgParts{p}, maskParts{p}, finalImgSize);
        
%         if maxImgSz > 256
%             imgParts{p} = imresize( imgParts{p}, 0.25);
%             maskParts{p} = imresize( maskParts{p}, 0.25);
%         end

        % do not take empty images 
        % (only background in mask or not wanted categories of events)
        categoriesInMaskPart = unique(maskParts{p});
        % take only in account categories of interest
        categoriesOfInterestInMaskPart = ...
            categoriesInMaskPart( ~ismember(categoriesInMaskPart, removeCategories) );
        
        if any(ismember(categoriesInMaskPart,removeCategories))
            keyboard
            figure
            imagesc(maskParts{p})
        end
             
        if numel(categoriesOfInterestInMaskPart) < 2
            continue
        end
        
        if sum(sum(maskParts{p}>0)) < 50
            continue
%             figure
%             imagesc(maskParts{p})
%             
%             keyboard
        end
        
        % new name for part of image
        newImgName = sprintf('imgData_%d.mat',imgNum);
        
        % set up path of created folder
        pathNewImg = fullfile(folderPath_splittedImgs, newImgName);

        % save image and mask as structure
        splittedImgData = struct( ...
            'img', imgParts{p}, ...
            'maskOfEvents', maskParts{p}, ...
            'pxSzX', pxSzX, ...
            'pxSzT', pxSzT );
        save(pathNewImg, 'splittedImgData')
        
        % save paths of splitted images
        splittedImgsPaths = [splittedImgsPaths; pathNewImg];
        
        % update img number
        imgNum = imgNum + 1;
    end
 
end

function [img_r, mask_r] = resizeImgToSpecificSize(img, mask, finalSz)

    % resize to chosen size
    img_r = imresize( img, finalSz/max(size(img)), 'Method', 'bilinear');
    mask_r = imresize( mask, finalSz/max(size(img)), 'Method', 'nearest');
    
end

















