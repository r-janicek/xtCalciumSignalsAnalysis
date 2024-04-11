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
stimFilePath_V = fullfile(stimulusFileList_voltage_subset.folder{ind}, ...
    stimulusFileList_voltage_subset.name{ind});
stimFilePath_caff = fullfile(stimulusFileList_voltage_subset.folder{ind}, ...
    strrep(stimulusFileList_voltage_subset.name{ind}, '_U_', '_S_'));
% delay to trigger confocal microscope
trigg_delay = seconds(str2double(hObjs.h_edit_triggerDelay.String));
% get igor file data
% stimulation protocol
stimulusData_V = IBWread(stimFilePath_V);
stimulus_img_t = seconds(linspace(stimulusData_V.x0, stimulusData_V.x1, ...
    stimulusData_V.Nsam)./1000) + stimulusFileList_voltage_subset.date(ind) - ...
    seconds((stimulusData_V.x1-stimulusData_V.x0)/1000);
% caffeine perfusion
stimulusData_caff = IBWread(stimFilePath_caff);
if isempty(stimulusData_caff.y)
    stimulusData_caff.y = nan(size(stimulusData_V.y));
end
% smooth data
s_const_caff = round(range(stimulusData_caff.y,'all')/5);
if s_const_caff==0, s_const_caff = 1; end
stimulusData_caff.y = round(stimulusData_caff.y./s_const_caff, 1, "significant").*s_const_caff;
s_const_V = round(range(stimulusData_V.y,'all')/5);
if s_const_V==0, s_const_V = 1; end
stimulusData_V.y = round(stimulusData_V.y./s_const_V, 1, "significant").*s_const_V;
if range(stimulusData_V.y)<5
    stimulusData_V.y = nan(size(stimulusData_V.y));
end

% start and end of image recording
imgRecStart = stimulusFileList_voltage_subset.date(ind) - ...
    seconds((stimulusData_V.x1-stimulusData_V.x0)/1000) + trigg_delay;
imgRecEnd = imgRecStart + ...
    seconds(imgData.pxSzT*size(imgData.wholeImgFluoXT,2)/1000);
% get stimulus during image recording
m_stimulus_img = stimulus_img_t<imgRecEnd & stimulus_img_t>imgRecStart;
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
figure('Position',[5 1*scSz(4)/4 scSz(3)-10 3*scSz(4)/4])
tiledlayout(4,1,"TileSpacing","tight")
% plot stimulation protocol
nexttile
plot(stimulus_img_t(m_stimulus_img), ...
    stimulusData_V.y(m_stimulus_img), 'LineWidth',3)
set(gca, 'FontSize',22, ...
    'XLim', [tempData_img.time(1) tempData_img.time(end)], 'XTick',[])
ylabel({"field stimulation","protocol"}, "FontSize",24)
% plot caffeine local perfusion
nexttile
plot(stimulus_img_t(m_stimulus_img), ...
    stimulusData_caff.y(m_stimulus_img), 'LineWidth',3)
set(gca, 'FontSize',22, ...
    'XLim', [tempData_img.time(1) tempData_img.time(end)], 'XTick',[])
ylabel("caffeine perfusion", "FontSize",24)
% plot temperature
nexttile
plot(tempData_img.time, tempData_img.temperature, 'LineWidth',3)
set(gca, 'FontSize',22, ...
    'XLim', [tempData_img.time(1) tempData_img.time(end)], ...
    'YLim', [floor(min(tempData_img.temperature(:))) ...
             ceil(max(tempData_img.temperature(:)))], ...
    'XTick',[])
ylabel("temperature (\circC)", "FontSize",24)
% show image
nexttile
imagesc(imgData.wholeImgFluoXT, ...
    'YData',[1 size(imgData.wholeImgFluoXT,1)], ...
    'XData',[tempData_img.time(1) tempData_img.time(end)])
set(gca, 'FontSize',22, 'YTick',[], ...
    'XLim', [tempData_img.time(1) tempData_img.time(end)])
xlabel("time")

set(mainFig, "Pointer","arrow")

end