function splittedImgsPaths = ...
    splitImagesInToParts(folderPath,imgDataPaths)    % (folderPath,imgDataPaths)

% create new folder with cutted/splitted images
folderPath_splittedImgs = fullfile(folderPath,'splittedImgs');
if exist(folderPath_splittedImgs, 'dir')
    rmdir(folderPath_splittedImgs, 's')
end
mkdir(folderPath_splittedImgs)

% cut images so max lenght is 512 px
maxImgSz = 512;

% add waitbar 
hw = waitbar(0,'0 %; 0 s','Name','splitting images into parts ...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(hw,'canceling',0)


splittedImgsPaths = {};
imgNum = 1;
tIter = nan(numel(imgDataPaths),1);
for i = 1:numel(imgDataPaths)
    % get time per iteration
    tStartIter = tic;
    
    % load image
    data = bfopen(imgDataPaths{i});
    imgDataFluo = data{1,1}{1,1};
    % get meta data, all of them
    omeMeta = data{1,4};
    % meta data
    metaData = cat(2,cell(data{1,2}.keySet.toArray),cell(data{1,2}.values.toArray));
    indx_SM = find(~cellfun(@isempty,regexp(metaData(:,1),'\<ScanMod+e','match'))==1,1,'first');
    for m = 1:size(metaData,1)
        ImageDescription_Names{m,1} = genvarname(char(metaData{m,1}));
        ImageDescription_val {m,1} = metaData{m,2};
    end
    [~,ia,~] = unique(ImageDescription_Names);
    metaDataS = cell2struct(cell(ImageDescription_val(ia)), cell(ImageDescription_Names(ia)),1);
    %
    pxSzT = str2double(metaDataS.GlobalTimePerLine)/1000; % in ms
    pxSzX = omeMeta.getPixelsPhysicalSizeX(0).value;
    
    pixelSize = ome.units.quantity.Length( pxSzX, ...
        ome.units.UNITS.MICROMETER );
    pixelSizeT = ome.units.quantity.Time( java.lang.Double(pxSzT), ...
        ome.units.UNITS.MILLISECOND );
    
    % check size of image, if necessary to split it
    if any( size(imgDataFluo) > 1024 )
        % split image into parts
        a  = size(imgDataFluo, 1);
        b  = size(imgDataFluo, 2);
        
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
        
        imgParts = mat2cell(imgDataFluo, a_sz, b_sz);
        
    else
        imgParts = {imgDataFluo};
    end
    
    % do not save very small parts
    m_sz_imgParts = cellfun(@(x) min(size(x)),imgParts) < 100;
    imgParts = imgParts(~m_sz_imgParts);
        
    for p = 1:numel(imgParts)
        
        % new name for part of image
        newImgName = sprintf('img_%d.ome.tiff',imgNum);
        
        % set up path of created folder
        pathNewImg = fullfile(folderPath_splittedImgs, newImgName);
        
        % create new metadata
        plane = imgParts{p};
        metadata = createMinimalOMEXMLMetadata(plane);
        metadata.setPixelsPhysicalSizeX(pixelSize, 0);
        metadata.setPixelsPhysicalSizeY(pixelSize, 0);
        metadata.setPlaneDeltaT(pixelSizeT,0,0);
        
        % save image
        bfsave(plane, pathNewImg, 'metadata', metadata);
        
        % save paths of splitted images
        splittedImgsPaths = [splittedImgsPaths; pathNewImg];
        
        % update img number
        imgNum = imgNum + 1;
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
    
end

% delete waitbar
delete(hw), clearvars hw

end

