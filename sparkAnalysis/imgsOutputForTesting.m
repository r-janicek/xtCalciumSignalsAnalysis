function imgsOutputForTesting(hO,~,mainFig)
% imgs output to test with CaSparks software
% save multiple images from loaded one in main window, and txt file with
% metadata (pixel size)
img_ext = 'bmp';
maxImgSzT = 1000; % pixels
% get data
imgData = getappdata(mainFig, 'imgData');
hObs = getappdata(mainFig, 'hObjs');
if isempty(imgData)
    return
end
% ask where to save splitted images
% previous path
try
    prevPath = getappdata(mainFig, 'allSplittedImgsDir');
catch
end
if isempty(prevPath)
    prevPath = imgData.filePath;
end
allSplittedImgsDir = uigetdir(prevPath, ...
    'select folder to save splitted images');
if allSplittedImgsDir == 0
    return
end
setappdata(mainFig, 'allSplittedImgsDir', allSplittedImgsDir)
[~,img_n,~] = fileparts(imgData.fileName);
% create folder with image name
if exist(fullfile(allSplittedImgsDir,img_n), 'dir')
    listing = struct2table(dir(allSplittedImgsDir));
    dirNames = listing.name(listing.isdir);
    idx = sum(cellfun(@(x) contains(x,img_n), dirNames));
    splittedImgDir = fullfile( ...
        allSplittedImgsDir, sprintf('%s_%d', img_n, idx));
    mkdir(splittedImgDir)
else
    splittedImgDir = fullfile(allSplittedImgsDir, img_n);
    mkdir(splittedImgDir)
end
% split image if its size in time is more than 1000 pixels and save them in
% selected dir
% first normalize image in range (0,1) then rotate it and convert to 8bit
wholeImg = imrotate( ...
    im2uint8( ...
    rescale( ...
    imgData.imgDataXTfluoR-imgData.blank) ), -90);

if size(wholeImg,1) > maxImgSzT
    n_imgs = ceil(size(wholeImg,1)/maxImgSzT);
    for i = 1:n_imgs
        s = (i-1)*maxImgSzT+1;
        e = min(i*maxImgSzT,size(wholeImg,1));
        splittedImg = wholeImg(s:e, :);
        % save it in selected dir in selected format
        imwrite(splittedImg, ...
            fullfile(splittedImgDir, ...
                     sprintf('%s_p%d.%s',img_n,i,img_ext)), img_ext)
    end
else
    % save it in selected dir in selected format
    imwrite(wholeImg, ...
        fullfile(splittedImgDir, sprintf('%s.%s',img_n,img_ext), img_ext))
end
% save csv file with pixel size
pxSz = [{'x (µm)', imgData.pxSzX}; ...
              {'t (ms)', imgData.pxSzT}];
writecell(pxSz, ...
    fullfile(splittedImgDir,sprintf('%s.csv',img_n)))









