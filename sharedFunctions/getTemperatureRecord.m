function getTemperatureRecord(hO,E,mainFig)
%{
load directory with both temperature and stumulus files

%}
% set pointer
set(mainFig, "Pointer","watch")

% data 
hObjs = getappdata(mainFig, 'hObjs');
imgData = getappdata(mainFig, 'imgData');
% get temperature and stimulus files from folder
% get folder path to temperature and stimulation data
allFilesDir = uigetdir(imgData.filePath);
allFilesList = struct2table(dir(fullfile(allFilesDir,'**/*.*')));
allFilesList(allFilesList.isdir, :) = [];
% get temperature files, sort them and read all of them
tempFileList = allFilesList(contains(allFilesList.name, '.txt'), :);
[~,fn_temp,~] = fileparts(tempFileList.name);
tempFileList.tempStarRecTime = ...
    datetime(fn_temp, 'InputFormat','yyyyMMdd-HH-mm-ss');
tempFileList = sortrows(tempFileList, "tempStarRecTime", "ascend");
tempDataAll = table();
for i = 1:height(tempFileList)
    % read temperature file
    f_temp = fullfile(tempFileList.folder{i}, tempFileList.name{i});
    opts = detectImportOptions(f_temp);
    tempData = readtable(f_temp, opts);
    tempData.Properties.VariableNames = {'time', 'temperature'};
    % create final table
    tempDataAll = [tempDataAll ;tempData];
end
tempDataAll.Properties.VariableUnits = {'', 'celsius'};
tempDataAll.time = datetime(tempDataAll.time, ...
    'InputFormat','yyyyMMdd-HH-mm-ss');

% get stimulus files
stimulusFileList = allFilesList(contains(allFilesList.name, '.bwav'), :);
stimulusFileList.date = datetime(stimulusFileList{:,'date'}, ...
    'InputFormat','dd-MMM-yyyy HH:mm:ss');
stimulusFileList_voltage = stimulusFileList( ...
    contains(stimulusFileList.name, '_U_'), :);
% find the stimulus file which most likely was used to
% trigger line-scan image recording
stimulusFileList_voltage_subset = stimulusFileList_voltage( ...
    stimulusFileList_voltage.date > imgData.imgCaptureDate, :);
[~, ind] = min(abs(stimulusFileList_voltage_subset.date - imgData.imgCaptureDate));
stimFilePath = fullfile(stimulusFileList_voltage_subset.folder{ind}, ...
    stimulusFileList_voltage_subset.name{ind});
% delay to trigger confocal microscope
trigg_delay = seconds(str2double(hObjs.h_edit_triggerDelay.String));
% get igor file data
stimulusData = IBWread(stimFilePath);
stimulus_img_t = linspace(stimulusData.x0, stimulusData.x1, stimulusData.Nsam);
stimulus_img_V = stimulusData.y;

img_rec_s = round(trigg_delay/seconds(stimulusData.dx/1000));
img_rec_e = round(img_rec_s + ...
    imgData.pxSzT*size(imgData.wholeImgFluoXT,2)/stimulusData.dx);
if img_rec_e > stimulusData.Nsam
    img_rec_e = stimulusData.Nsam;
end

% figure
% plot(stimulus_img_V(img_rec_s:img_rec_e))

% start and end of image recording
imgRecStart = stimulusFileList_voltage_subset.date(ind) - ...
    seconds((stimulusData.x1-stimulusData.x0)/1000) + trigg_delay;
imgRecEnd = imgRecStart + ...
    seconds(imgData.pxSzT*size(imgData.wholeImgFluoXT,2)/1000);
% get temperature during image recording
m_temperature = tempDataAll.time<imgRecEnd & tempDataAll.time>imgRecStart;
tempData_img = tempDataAll(m_temperature, :);
[~, TFrm] = rmoutliers(tempData_img.temperature, ...
    "gesd");%, seconds(5), "SamplePoints",t);
tempData_img.temperature(TFrm) = nan;
tempData_img.temperature = fillmissing(tempData_img.temperature,'linear');

% save temperature data
imgData.tempData_img = tempData_img;
setappdata(mainFig, 'imgData', imgData)

% show figure with temperature recording and line scan image
scSz = get(0,'ScreenSize');
figure('Position',[5 3*scSz(4)/4 scSz(3)-10 scSz(4)/2])
tiledlayout(2,1,"TileSpacing","tight")
% plot temperature
nexttile
plot(tempData_img.time, tempData_img.temperature, 'LineWidth',3)
set(gca, 'FontSize',22, ...
    'XLim', [tempData_img.time(1) tempData_img.time(end)], ...
    'YLim', [floor(min(tempData_img.temperature(:))) ...
             ceil(max(tempData_img.temperature(:)))], ...
    'XTick',[])
ylabel("temperature (\circC)", "FontSize",24)
nexttile
imagesc(imgData.wholeImgFluoXT, ...
    'YData',[1 size(imgData.wholeImgFluoXT,1)], ...
    'XData',[tempData_img.time(1) tempData_img.time(end)])
set(gca, 'FontSize',22, 'YTick',[], ...
    'XLim', [tempData_img.time(1) tempData_img.time(end)])
xlabel("time")

set(mainFig, "Pointer","arrow")

end