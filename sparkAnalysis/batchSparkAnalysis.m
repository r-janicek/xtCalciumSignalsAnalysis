function batchSparkAnalysis(~,~,mainFig)

% handles to objects in main window
hObjs = getappdata(mainFig,'hObjs');

% set batch processing mode
setappdata(mainFig,'batchProcessing',true)

datasetForMaskRCNN = false;


%% get paths to image files
% load whole selected directory
folderPath = uigetdir( ...
    '/Users/rjanicek/Desktop/data/DKI_S2808A14A_project', ...
    'Analyze all images from the folder');
if folderPath == 0
    return
end

fileList = struct2table(dir(fullfile(folderPath,'**/*.*')));
fileList = fullfile(fileList.folder,fileList.name);
% take only images (with olympus extensions: .oif, .oib)
%imgDataPaths = fileList(~cellfun(@isempty,regexp(fileList,'(.oif(?!\.files)|.oib)')));
imgDataPaths = fileList(~cellfun(@isempty,regexp(fileList,'.bmp')));

% remove caffeine images
imgDataPaths = imgDataPaths(cellfun(@isempty,regexp(imgDataPaths,'caf')));


%% initialize stuff for creating train and val datasets for mask-RCNN 
if datasetForMaskRCNN
  
    imgDataPaths = splitImagesInToParts(folderPath,imgDataPaths);
    
    % info dataset
    infoDataset = struct( ...
        'description', 'calciumSignals', ...
        'url', '', ...
        'version', 1.0, ...
        'year', '2021', ...
        'contributor', 'Rado Janicek', ...
        'date_created', datestr(today('datetime')) );
    
    % licenses for dataset
    licenseDataset = table( ...
        {''}, 1, {'license'} ,...
        'VariableNames',{'url', 'id', 'name'});
    
    % categories
    categoriesInDataset = table( ...
        {'calciumSignals'}, 1, {'spark'} ,...
        'VariableNames',{'supercategory', 'id', 'name'});
  
%     table( ...
%         {'calciumSignals';'calciumSignals';'calciumSignals';'calciumSignals'},...
%         [1;2;3;4], ...
%         {'spark'; 'longLastingSpark'; 'miniWave'; 'wave'} ,...
%         'VariableNames',{'supercategory', 'id', 'name'})

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
    
    % create folders for train and val
    folderPath_train = fullfile(folderPath,'train');
    folderPath_val = fullfile(folderPath,'val');
    
    if exist(folderPath_train, 'dir')
        rmdir(folderPath_train, 's')
    end
    if exist(folderPath_val, 'dir')    
        rmdir(folderPath_val, 's')
    end
    
    mkdir(folderPath_train)
    mkdir(folderPath_val)
    
    % annotation id, starts at 10000
    annotID = 10000;
end


%% main loop, get description from paths of images and analyze images for sparks
% change pointer to watch
set(mainFig,'Pointer','watch')
drawnow
% add waitbar 
hw = waitbar(0,'0 %; 0 s','Name','analyzing data ...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(hw,'canceling',0)

notAnalyzed = {};
allSparks = table();
tIter = nan(numel(imgDataPaths),1);
for i = 1:numel(imgDataPaths)
    try
        % get time per iteration
        tStartIter = tic;
        
        % extract description of experiment from path of image
        expInfo = getExperimentInfo(imgDataPaths{i});
        
        % get animal and pharmacology from excel file
        % specific for data for IP3 project
        projectName = "IP3 local signals";
        projectName = "sparks PA-RFP-RyR";
        expInfo.animal = {'PA-RFP RyR'};
        if contains(expInfo.imgName,'ctrl')
            expInfo.expCond = {'ctrl'};
        elseif contains(expInfo.imgName,'iso')
            expInfo.expCond = {'iso'};
        else
            expInfo.expCond = {'NA'};
        end
         

        
%         if strcmp(expInfo.expCond,'blank only')
%             continue
%         else
%             % get animal and pharmacology from excel file
%             expCondFile = ...
%                 fileList{contains(fileList, 'expConditionsAllData.xlsx')};
%             expCondData_cell = readcell(expCondFile);
%             expCondData = cell2table(expCondData_cell(2:end,1:5), ...
%                 'VariableNames', expCondData_cell(1,1:5));
%             % get experimental conditions for individual images
%             m_expCond = expCondData.date == str2double(expInfo.date) & ...
%                 strcmp(expCondData.cell, lower(expInfo.cell)) & ...
%                 expCondData.protocolNum == str2double(expInfo.t_expCond);
%             
%             expInfo.animal = expCondData.animal(m_expCond);
%             expInfo.expCond = expCondData.condition(m_expCond);
%             % get blank from image number 1 from the cell
%             blankImg_path = fileList( contains(fileList,expInfo.date) & ...
%                 contains(fileList, ['_',expInfo.cell{1},'1','_']));
%             blankImgData = bfopen(blankImg_path{1});
%             blankImg = blankImgData{1,1}{1,1};
%             % treshlod image and get blank from out of cell area
%             T = graythresh(blankImg);
%             BW = imbinarize(blankImg,T);
%             se = strel('disk',5);
%             BW = imerode(imclose(BW,se),se);
%             % try to find better method to get BW mask
% %             figure
% %             imshowpair(blankImg, BW, 'montage')
%             blank = round(mean(blankImg(~BW),"all"),3);            
%         end
     
        % load image
        setappdata(mainFig,'FilePath',expInfo.imgPath)
        setappdata(mainFig,'FileName',expInfo.imgName)
        openImg([],[],mainFig)

        % get image data
        imgData = getappdata(mainFig,'imgData');
        % set up blank
        experimentalApproach = 'permeabilized myocytes';
        if strcmp(experimentalApproach, 'permeabilized myocytes')
            if ~isempty(imgData.userName)
                switch imgData.userName
                    case 'Rado'
                        blank = 190;
                    case 'Duilio'
                        blank = 152;
                end
            else
               blank = 0;
            end
        end
        if exist("blank","var")
            hObjs.h_edit_Blank.String = num2str(blank);
            % estimation of background in image in analysis and crop out
            % background
            % calculate x-profile
            % x_prof = mean(imgData.imgDataXTfluoF,2);
            % t_prof = mean(imgData.imgDataXTfluoF,1);
            % cellDetThrs = blank*1.25;
            % crop_x1 = find(x_prof>cellDetThrs,1,'first');
            % crop_x2 = numel(x_prof) - ...
            %     find(flipud(x_prof)>cellDetThrs,1,'first') + 1;
            % crop_t1 = find(t_prof>cellDetThrs,1,'first');
            % crop_t2 = numel(t_prof) - ...
            %     find(flipud(t_prof)>cellDetThrs,1,'first') + 1;
            % % try to look for flash position
            % if contains(expInfo.expCond, 'UVflash', 'IgnoreCase',true)
            %     if ~isempty(imgData.imgDataXTtrans)
            %         [~,flash_pos] = max( ...
            %             gradient( ...
            %             smooth(mean(imgData.imgDataXTtrans,1), 5) ) );
            %     else
            %         [~,flash_pos] = max(gradient(smooth(t_prof,5)));
            %     end
            %     flash_pos = flash_pos + 10;
            %     if isempty(flash_pos)
            %         flash_pos = 1;
            %     end
            %     if flash_pos<crop_t1
            %         flash_pos = crop_t1;
            %     end
            %     if flash_pos>crop_t2
            %         flash_pos = crop_t1;
            %     end
            %     crop_t1 = flash_pos;
            % 
            % end

            % call crop function
            % cropImg(hObjs.h_pb_crop,[],mainFig)
            % % set up rectangle to crop image
            % h_rect_crop = getappdata(mainFig,'h_rect_crop');
            % % new pos of rectangle
            % h_rect_crop_newPos(1) = crop_t1*imgData.pxSzT;
            % h_rect_crop_newPos(2) = crop_x1;
            % h_rect_crop_newPos(3) = (crop_t2-crop_t1)*imgData.pxSzT;
            % h_rect_crop_newPos(4) = crop_x2-crop_x1;
            % h_rect_crop.setPosition(h_rect_crop_newPos)
            % % crop image
            % cropImg(hObjs.h_pb_crop,[],mainFig)

        else
            % check if there are no zeros in image, if so than crop it
            % get 10th percentile from xy image - estimation of background
            if ~isempty(imgData.imgDataXY)
                bg_est = double(prctile(imgData.imgDataXY{1}(:),5));
                if prctile(imgData.wholeImgFluoXT(:),10) < bg_est*1.05

                    % add better check for line between cell and background
                    mImg = imgData.imgDataXTfluoFN - bg_est*1.01;
                    mImg(mImg>0) = 1;
                    mImg(mImg<=0) = 0;

                    % find biggest rectangle, which is in image area
                    img_CC_event = bwconncomp(logical(mImg),8);

                    %calculate properties of all CC (regions)
                    statEvents = regionprops(img_CC_event,imgData.imgDataXTfluoFN, ...
                        'SubarrayIdx', 'Area','BoundingBox');
                    [~,indMaxArea] = max([statEvents.Area]);

                    % call crop function
                    cropImg(hObjs.h_pb_crop,[],mainFig)
                    % set up rectangle to crop image
                    h_rect_crop = getappdata(mainFig,'h_rect_crop');
                    % new pos of rectangle
                    h_rect_crop_newPos = statEvents(indMaxArea).BoundingBox;
                    h_rect_crop_newPos(1) = h_rect_crop_newPos(1)*imgData.pxSzT;
                    h_rect_crop_newPos(3) = h_rect_crop_newPos(3)*imgData.pxSzT;
                    h_rect_crop.setPosition(h_rect_crop_newPos)
                    % crop image
                    cropImg(hObjs.h_pb_crop,[],mainFig)
                end
            end
        end
        % set up main window
        switch expInfo.animal{1}
            case 'DKI'
                animalStrain = 'S2808A/14A';
                
            case 'S2808A14A'
                animalStrain = 'S2808A/14A';
                
            case 'WT'
                animalStrain = 'wt';

            case 'IP3_OE'
                animalStrain = 'IP3_OE';

            case {'IP3R/tTA', 'IP3R:tTA'}
                animalStrain = 'IP3R/tTA';
                
            otherwise
                animalStrain = expInfo.animal{1};
        end
        
        % set up animal
        hObjs.popUpMenuAnimal.Value = ...
            find(strcmpi(hObjs.popUpMenuAnimal.String,animalStrain));
    
        % set up notes
        expNotes = hObjs.h_table_notes.Data;
        expNotes(1) = expInfo.expCond;
        if strcmp(projectName, "IP3")
            expNotes(2:end) = [];
        elseif strcmp(projectName, "IP3 local signals")

        elseif strcmp(projectName, "sparks PA-RFP-RyR")
            expNotes(1) = expInfo.expCond;
            hObjs.h_edit_pxSzT.String = num2str(2.1);
            hObjs.h_edit_pxSzX.String = num2str(0.138);
            setPxSzManually(hObjs.h_edit_pxSzT, [], mainFig)
            setPxSzManually(hObjs.h_edit_pxSzX, [], mainFig)
        else
            expNotes(4) = {sprintf('apl. of drug t = %d (minutes)',expInfo.t_expCond)};
        end
        set(hObjs.h_table_notes,'Data',expNotes)
        
        % find sparks
        eventsDetection([], [], mainFig, 'new')
        
        if datasetForMaskRCNN
            % read imgData, with updated structure
            imgData = getappdata(mainFig,'imgData');
            caEvents = getappdata(mainFig,'sparkDetection');
            
            % remove sparks with certain properties
            maskDetectedEvents = true(1,numel(caEvents.detectedEvents));
            for e = 1:numel(caEvents.detectedEvents)
                imgEvent = imgData.imgDataXTfluoFN( ...
                    caEvents.detectedEvents(e).SubarrayIdx{1}, ...
                    caEvents.detectedEvents(e).SubarrayIdx{2} );
                % filter
               
                imgEvent = imgaussfilt(imgEvent,1);
                % check if 5% of pixels of events are higher than 2*noise
                if (prctile(imgEvent(caEvents.detectedEvents(e).Image), 90) < ...
                        2*imgData.stdNoise) || ...
                     size(imgEvent,2)* imgData.pxSzT < str2double(hObjs.h_edit_MinDurSpark.String) || ...
                     size(imgEvent,1)* imgData.pxSzX < str2double(hObjs.h_edit_MinWidthSpark.String)
                    
                   maskDetectedEvents(e) = false; 
                   caEvents.detectedEventsRec(e).EdgeColor = [0 0 0]; 
                end
 
            end
            
            detectedEventsFiltered = ...
                caEvents.detectedEvents(maskDetectedEvents);
            
            if isempty(detectedEventsFiltered)
                continue
            end
            
            % image structure
            img = struct( ...
                'id', i, ...
                'license', licenseDataset.id, ...
                'coco_url', '', ... 
                'flickr_url', '', ...
                'width', size(imgData.imgDataXTfluoFN,2), ...
                'height', size(imgData.imgDataXTfluoFN,1), ...
                'file_name', sprintf('img_%d.png',i), ...
                'date_captured', '');
            
            % create structure for all anotations in image
            allAnnotationsInImg = struct( ...
                'id', {}, ...
                'category_id', {}, ...
                'iscrowd', {}, ...
                'segmentation', {}, ... % [x1,y1,x2,y2,...]
                'image_id', {}, ...
                'area', {}, ... % number of pixels
                'bbox', {} );  %  boundig box [x,y,w,h] [0 0] is left upper corner
            for j = 1:numel(detectedEventsFiltered)
         
                eventBorder = detectedEventsFiltered(j).ConvexHull;
                bbox = detectedEventsFiltered(j).BoundingBox;
                % in python idx starts at 0
                bbox = [ bbox(1)-1, bbox(2)-1, bbox(3), bbox(4) ]; 
                eventBorder = eventBorder-1;
                eventBorder = eventBorder';
                eventBorder = eventBorder(:);
                eventBorder = eventBorder';
                
                singleAnnotationsInImg = struct( ...
                    'id', annotID, ...
                    'category_id', 1, ...
                    'iscrowd', 0, ...
                    'segmentation', {{eventBorder}}, ... % [x1,y1,x2,y2,...]
                    'image_id', i, ...
                    'area', detectedEventsFiltered(j).Area, ... % number of pixels
                    'bbox', bbox );  %  boundig box [x,y,w,h] [0 0] is left upper corner

                allAnnotationsInImg = ...
                    [allAnnotationsInImg; singleAnnotationsInImg];
                % add 1 to annotaion ID
                annotID = annotID+1;
            end
            
            
            % first 80% of images as train dataset, 20% as validation
            % dataset
            if i <= round(numel(imgDataPaths)*0.8)
            	% save image as jpeg
                imwrite(uint16(imgData.imgDataXTfluoR), ...
                    fullfile(folderPath_train,img.file_name), ...
                    'BitDepth',16);
                
                trainImgs = [trainImgs; img];
                annotationsInTrainImgs = ...
                    [annotationsInTrainImgs; allAnnotationsInImg];
                
            else
                % save image as jpeg
                imwrite(uint16(imgData.imgDataXTfluoR), ...
                    fullfile(folderPath_val,img.file_name), ...
                    'BitDepth',16);
                
                valImgs = [valImgs ; img];
                annotationsInValImgs = ...
                    [annotationsInValImgs; allAnnotationsInImg];
                
            end
            
            
        else
            % filter and analyze sparks
            detectedSparksAnalysis([],[],mainFig)
      
            % get analyzed data
            sparkDetection = getappdata(mainFig,'sparkDetection');
            imgData = getappdata(mainFig,'imgData');
            imgSzX = size(imgData.imgDataXTfluoFN,1)*imgData.pxSzX;        % um
            imgSzT = size(imgData.imgDataXTfluoFN,2)*imgData.pxSzT-imgData.pxSzT;  % ms
            
            % create output allSparks table
            dataImgSparks = struct2table(sparkDetection.eventParams);
            dataImgSparks.maskOfAcceptedSparks = sparkDetection.maskOfAcceptedSparks(:);
            dataImgSparks.imgSzX = repmat( imgSzX, [height(dataImgSparks),1] );
            dataImgSparks.imgSzT = repmat( imgSzT, [height(dataImgSparks),1] );
            
            imgInfo = table( ...
                repmat( {expInfo.imgPath}, [height(dataImgSparks),1] ), ...
                repmat( expInfo.animal, [height(dataImgSparks),1] ), ...
                repmat( expInfo.date, [height(dataImgSparks),1] ), ...
                repmat( expInfo.loadingGroup, [height(dataImgSparks),1] ), ...
                repmat( expInfo.cell, [height(dataImgSparks),1] ), ...
                repmat( expInfo.expCond, [height(dataImgSparks),1] ), ...
                repmat( expInfo.imgName, [height(dataImgSparks),1] ), ...
                repmat( expInfo.t_expCond, [height(dataImgSparks),1] ), ...
                repmat( sparkDetection.sparkFreq, [height(dataImgSparks),1] ), ...
                repmat( sparkDetection.correctedSparkFreq, [height(dataImgSparks),1] ), ...
                'VariableNames', {'pathImg' 'animal' 'day' 'loading' 'cell' ...
                'expCond' 'imgName' 'drugAplTime' 'sparkFreq' 'correctedSparkFreq'});
            
            allSparks = [ allSparks; [imgInfo,dataImgSparks] ];
            % save results
            % saveSparkDetAnalysis([],[],mainFig)
            
        end
        
        % get time per iteration, in seconds
        tIter(i) = toc(tStartIter);
        estimTimeOfFinish = mean(tIter,'omitnan')*(numel(imgDataPaths)-i);
        % in hours
        estimTimeOfFinish_hour = [floor(estimTimeOfFinish/3600), ...
            round(mod(estimTimeOfFinish,3600)/60)]
        % in minutes
        estimTimeOfFinish_min = floor(estimTimeOfFinish/60)
        % in seconds
        estimTimeOfFinish_sec = round(estimTimeOfFinish) % round(mod(estimTimeOfFinish,60));
        
        if estimTimeOfFinish_hour(1) >= 1
            % Report current estimate in the waitbar's message field
            waitbar( i/numel(imgDataPaths),hw,...
                sprintf( '%d %%; remaining time: %d hours and %d minutes', ...
                round(100*i/numel(imgDataPaths)),...
                estimTimeOfFinish_hour(1),estimTimeOfFinish_hour(2) ) )
        else
            if estimTimeOfFinish_min > 1
                % Report current estimate in the waitbar's message field
                waitbar( i/numel(imgDataPaths),hw,...
                    sprintf( '%d %%; remaining time: %d minutes', ...
                    round(100*i/numel(imgDataPaths)),...
                    estimTimeOfFinish_min ) )
            else
                % Report current estimate in the waitbar's message field
                waitbar( i/numel(imgDataPaths),hw,...
                    sprintf( '%d %%; remaining time: %d seconds', ...
                    round(100*i/numel(imgDataPaths)),...
                    estimTimeOfFinish_sec ) )
            end
        end
        
        % Check for clicked Cancel button
        if getappdata(hw,'canceling')
            % delete waitbar
            delete(hw)
            break
        end
        
    catch
        notAnalyzed = [notAnalyzed;imgDataPaths{i}];
    end
    
end

% delete waitbar
delete(hw), clearvars hw


if datasetForMaskRCNN
    % create and save final json file with annotations
    % with name of folder and save it into folder
    final_json_file_train = jsonencode( struct( ...
                'info', infoDataset, ...
                'licenses', licenseDataset, ...
                'images', trainImgs, ...
                'annotations', annotationsInTrainImgs, ...
                'categories', categoriesInDataset ) );
       
    final_json_file_val = jsonencode( struct( ...
                'info', infoDataset, ...
                'licenses', licenseDataset, ...
                'images', valImgs, ...
                'annotations', annotationsInValImgs, ...
                'categories', categoriesInDataset ) );

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
    
    % put back pointer to arrow
    set(mainFig,'Pointer','arrow')
    drawnow   
    return
    
else
    % set up units of table
    allSparks.Properties.VariableUnits = ...
        {'' '' '' '' '' '' '' 'minutes' ...
        '#sp*100um-1*s-1' '#sp*100um-1*s-1' '' ...
        'deltaF/F0' 'ms' 'ms' 'ms' 'ms' 'ms' 'um' 'deltaF/F0*um^3' ...
        'ms' 'ms' 'deltaF/F0*ms' '' '' '' '' '' '' 'um' 'ms'};


end
% put back pointer to arrow
set(mainFig,'Pointer','arrow')
drawnow        
 

%% save all sparks table 
% set up path for xls file
path_xls = fullfile(folderPath,'sparkDetResults_Matlab.xls');

% delete xls file if exist
if exist(path_xls, 'file')
    delete(path_xls)
end

% prepare data for writing into excel
analysisOf = {sprintf('analyzed sparks from folder: %s',folderPath)};

% parameters of analysis
parametersOfSparkDetection = ...
    [ [{'paramsOfSparkDetection:'},{''}];
    [{'fFWHM (px)'},str2double(hObjs.h_edit_fFWHM.String)];
    [{'smoothIter'},str2double(hObjs.h_edit_smoothIter.String)];
    [{'baseIter'},str2double(hObjs.h_edit_baseIter.String)];
    [{'tresh'},str2double(hObjs.h_edit_tresh.String)];
    [{'expFactor'},str2double(hObjs.h_edit_expFactor.String)];
    [{'minSparkDuration (ms)'},str2double(hObjs.h_edit_MinDurSpark.String)];
    [{'minSparkWidth (um)'},str2double(hObjs.h_edit_MinWidthSpark.String)];
    [hObjs.txt_spDet.String,{''}] ];

sparkSelectionCriteria = ...
    [ [{'method to analyze sparks parameters:'},...
         hObjs.calcSpParamMethodRBgroup.SelectedObject.String ];
    [{'sparkSelectionCriteria:'},{''} ];
    [{'spAmpl > # (F/F0)'},str2double(hObjs.h_edit_noise.String) ];
    [{'spFDHM > # (ms)'},str2double(hObjs.h_edit_spFDHM.String) ];
    [{'spFWHM > # (um)'},str2double(hObjs.h_edit_spFWHM.String) ];
    [{'smoothingSpan (ms)'},str2double(hObjs.h_smooth_edit.String) ];
    [{'baseline detection sensitivity'},str2double(hObjs.h_bsDet_edit.String) ] ];
    
% create sparks data table
resultSparksParamsNames = cellfun(@(x,y) [x,' (',y,')'], ...
    allSparks.Properties.VariableNames, ...
    allSparks.Properties.VariableUnits, ...
    'UniformOutput',0);

resultSparksParams = [resultSparksParamsNames;table2cell(allSparks)];

% save analysis info and detected sparks paramaters
infoCellArr = [ [analysisOf, cell(1,1)]; cell(1,2); ...
                parametersOfSparkDetection; cell(1,2); ...
                sparkSelectionCriteria];
writecell(infoCellArr,path_xls,'Sheet','info');

% check for the length of result table if it not too much for .xls (more than 60000)
% than split it
maxRowsOfData = 60000;
if size(resultSparksParams,1)-1 > maxRowsOfData
    % multiple sheets of results
    varNames = resultSparksParams(1,:);
    data = resultSparksParams(2:end,:);
    
    resSheetsNum = ceil((size(data,1)/maxRowsOfData));
    
    br = [(0:resSheetsNum-1)*maxRowsOfData+1, size(data,1)+1];
    for s = 1:resSheetsNum
        partOfDataToWrite = [ varNames; data(br(s):br(s+1)-1,:)];
        writecell(partOfDataToWrite,path_xls,'Sheet',sprintf('results_%d',s));
        %xlwrite(path_xls, partOfDataToWrite, sprintf('results_%d',s), 'A1');
    end
else
    % single sheet of results
    %xlwrite(path_xls, resultSparksParams, 'results', 'A1');
    writecell(resultSparksParams,path_xls,'Sheet','results');
end

end


%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
function expInfo = getExperimentInfo(imgDataPath)

% try to find animal and cell
infoStr = split(imgDataPath,'/');
[pImg,nImg,extImg] = fileparts(imgDataPath);

% check for time in file name, also test for if there is folder loading
try
    if contains(nImg,'Tamara')
        % specific for data from IP3 project
        % add path
        expInfo.imgPath = pImg;
        % get cell name and protocol number
        % imgIdx = regexpi(nImg,'(?<=[_,-])\D+\d+(?=_raw)','match');
        cellAndProt = regexpi(nImg, ...
            '(?<=[_,-])(?<cell>\D+)(?<protocolNum>\d+)(?=_raw)','names');
        if str2double(cellAndProt.protocolNum)==1
            % image of the cell, use only to get blank
            expInfo.expCond = {'blank only'};
        else
            expInfo.expCond = {''};
        end
        expInfo.animal = {'wt'};
        expInfo.date = regexpi(pImg, '(?<=/)\d+(?=/)', 'match');
        expInfo.loadingGroup = {'NA'};
        expInfo.cell = {cellAndProt.cell};
        expInfo.imgName = {[nImg,extImg]};
        expInfo.t_expCond = {cellAndProt.protocolNum};

    else
        % add path
        expInfo.imgPath = pImg;

        % possible animal strain names
        possibleAnimNames = {'wt','S2808','S2808A', ...
            'S2814','S2814A','S2030','S2030A',...
            'S2808A/14A','S2808A14A','R420Q','IP3R/tTA'};
        possibleAnimNames = strrep(possibleAnimNames, '/', ':');

        for animN = possibleAnimNames
            tf = strcmpi(infoStr,animN);
            if any(tf)
                indAnimInPath = find(tf,1,'last');
                expInfo.animal = animN;
                break
            end
        end

        % part of path from the string of animal
        partOfPath = infoStr(indAnimInPath:end);

        % add date
        expInfo.date = partOfPath(2);

        % add loading, if present
        indLoad = find( cellfun(@(x) ~isempty(regexp(x,'permeab|loading')), partOfPath), ...
            1,'first');
        if isempty(indLoad)
            expInfo.loadingGroup = {'NA'};
        else
            expInfo.loadingGroup = partOfPath(indLoad);
        end
        % check for cell string
        indCell = find( cellfun(@(x) ~isempty(regexp(x,'cel.\d*')), partOfPath), ...
            1,'first');
        if isempty(indCell)
            expInfo.cell =  'cell1';
            indCell = 1;
        else
            expInfo.cell =  partOfPath(indCell);
        end
        % expInfo.expCond =  partOfPath(indCell+1);
        % get condition from image name
        expCond = split(nImg, '_');
        expInfo.expCond =  join(expCond(2:end),'-');
        expInfo.imgName = {[nImg,extImg]};
        if strcmpi(partOfPath{indCell+1},'PP1')
            % get time after application of drug, if any
            %str_expCond = regexp(nImg,'[_-]\d+(m)?','match');
            nImg_NoDot = regexp(nImg,'\.','split');
            str_expCond = regexp(nImg_NoDot{1},'\d+','match');
            t_expCond = str2double(str_expCond{end});
            expInfo.t_expCond = t_expCond;
        else
            expInfo.t_expCond = nan;
        end

    end

catch
    expInfo.imgPath = pImg;
    expInfo.animal = {'wt'};
    expInfo.date = {'NA'};
    expInfo.loadingGroup = {'NA'};
    expInfo.cell = {'NA'};
    expInfo.expCond = {'NA'};
    expInfo.imgName = {[nImg,extImg]};
    expInfo.t_expCond = nan;
end

end


