function getTemperatureRecord(hO, E, mainFig, ...
    selectPathToDir, showLoadedData)
%{
load directory with both temperature and stimulus files
parameters of function:
hO = handle to object
E = event structure
mainFig = handle to main figure
showLoadedData
selectPathToDir
%}

% set pointer
set(mainFig, "Pointer","watch")

% data 
hObjs = getappdata(mainFig, 'hObjs');
imgData = getappdata(mainFig, 'imgData');
% get temperature and stimulus files from folder
% get folder path to temperature and stimulation data
if selectPathToDir
    allFilesDir = uigetdir(imgData.filePath);
    hObjs.h_text_tempAndStimDirPath.String = allFilesDir;
else
    allFilesDir = hObjs.h_text_tempAndStimDirPath.String;
end
allFilesList = struct2table(dir(fullfile(allFilesDir,'**/*.*')));
if isempty(allFilesList)
    opts = struct('WindowStyle','modal',... 
              'Interpreter','tex');
    errorDlg = errordlg( ...
        '\fontsize{20} Check folder path to load temperature and stimulus files!', ...
        'Folder Error', opts);
    figure(errorDlg)
    return
end
allFilesList(allFilesList.isdir, :) = [];
% get temperature files, sort them and read all of them
tempFileList = allFilesList(contains(allFilesList.name, '.txt'), :);
[~,fn_temp,~] = fileparts(tempFileList.name);
tempFileList.tempStarRecTime = ...
    datetime(fn_temp, 'InputFormat','yyyyMMdd-HH-mm-ss', ...
    'Format','yyyy-MM-dd HH:mm:ss.SSS');
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
    'InputFormat','yyyyMMdd-HH-mm-ss', ...
    'Format','yyyy-MM-dd HH:mm:ss.SSS');

% get stimulus files
stimulusFileList = allFilesList(contains(allFilesList.name, '.bwav'), :);
stimulusFileList.date = datetime(stimulusFileList{:,'date'}, ...
    'InputFormat','dd-MMM-yyyy HH:mm:ss', ...
    'Format','yyyy-MM-dd HH:mm:ss.SSS');
stimulusFileList_voltage = stimulusFileList( ...
    contains(stimulusFileList.name, '_U_'), :);
% find the stimulus file which most likely was used to
% trigger line-scan image recording
stimulusFileList_voltage_subset = stimulusFileList_voltage( ...
    stimulusFileList_voltage.date > imgData.imgCaptureDate, :);
[~, ind] = min(abs(stimulusFileList_voltage_subset.date - imgData.imgCaptureDate));
try 
    timeDiff = stimulusFileList_voltage_subset.date(ind) - ...
        imgData.imgCaptureDate;
    if timeDiff > minutes(5)
        falseFile = true;
    else
        falseFile = false;
    end
catch
    falseFile = true;
end
if isempty(ind) || falseFile
    opts = struct('WindowStyle','modal',...
        'Interpreter','tex');
    errorDlg = errordlg( ...
        '\fontsize{20} Check folder path to load temperature and stimulus files!', ...
        'Could not find proper stimulus file', opts);
    figure(errorDlg)
    return
end
stimFilePath_V = fullfile(stimulusFileList_voltage_subset.folder{ind}, ...
    stimulusFileList_voltage_subset.name{ind});
stimFilePath_caff = fullfile(stimulusFileList_voltage_subset.folder{ind}, ...
    strrep(stimulusFileList_voltage_subset.name{ind}, '_U_', '_S_'));
% delay to trigger confocal microscope
trigg_delay = seconds(str2double(hObjs.h_edit_triggerDelay.String));
% get igor file data
% stimulation protocol
stimulusData_V = IBWread(stimFilePath_V);
% downsample to 1 ms
n_per_ms = round(1/stimulusData_V.dx);
stimulusData_V.y = stimulusData_V.y(1:n_per_ms:end);
stimulusData_V.Nsam = numel(stimulusData_V.y);
% in ms
stimulus_img_t = linspace(stimulusData_V.x0, stimulusData_V.x1, ...
    stimulusData_V.Nsam);
% caffeine perfusion
stimulusData_caff = IBWread(stimFilePath_caff);
if isempty(stimulusData_caff.y)
    stimulusData_caff.y = nan(size(stimulusData_V.y));
end
% downsample to 1 ms
stimulusData_caff.y = stimulusData_caff.y(1:n_per_ms:end);
% smooth data
s_const_caff = round(range(stimulusData_caff.y,'all')/5);
if s_const_caff==0, s_const_caff = 1; end
stimulusData_caff.y = round(stimulusData_caff.y./s_const_caff, 1, "significant").*s_const_caff;
% round around 0
[~, hist_caff_edges] = histcounts(stimulusData_caff.y);
stimulusData_caff.y( ...
    stimulusData_caff.y >= hist_caff_edges(find(hist_caff_edges<=0, 1, 'first')) & ...
    stimulusData_caff.y <= hist_caff_edges(find(hist_caff_edges>0, 1, 'first'))) = 0;
% stimulus
% smooth data
s_const_V = round(range(stimulusData_V.y,'all')/5);
if s_const_V==0, s_const_V = 1; end
stimulusData_V.y = round(stimulusData_V.y./s_const_V, 1, "significant").*s_const_V;
% round around 0
[~, hist_V_edges] = histcounts(stimulusData_V.y);
stimulusData_V.y( ...
    stimulusData_V.y >= hist_V_edges(find(hist_V_edges<=0, 1, 'first')) & ...
    stimulusData_V.y <= hist_V_edges(find(hist_V_edges>0, 1, 'first'))) = 0;
% if range(stimulusData_V.y)<5
%     stimulusData_V.y = nan(size(stimulusData_V.y));
% end

% start and end of image recording
imgRecStart = stimulusFileList_voltage_subset.date(ind) - ...
    milliseconds((stimulusData_V.x1-stimulusData_V.x0)) + trigg_delay;
imgRecEnd = imgRecStart + ...
    milliseconds(imgData.pxSzT*size(imgData.wholeImgFluoXT,2));
% get temperature during image recording
m_temperature = tempDataAll.time<imgRecEnd & tempDataAll.time>imgRecStart;
tempData_img = tempDataAll(m_temperature, :);
[~, TFrm] = rmoutliers(tempData_img.temperature, ...
    "gesd");%, seconds(5), "SamplePoints",t);
tempData_img.temperature(TFrm) = nan;
tempData_img.temperature = fillmissing(tempData_img.temperature,'linear');
% check if there are some temperature data
if isempty(tempData_img)
    opts = struct('WindowStyle','modal',...
        'Interpreter','tex');
    errorDlg = errordlg( ...
        '\fontsize{20} Check folder path to load temperature and stimulus files!', ...
        'Could not find proper temperature file', opts);
    figure(errorDlg)
    return
end
% time of temperature recording in milliseconds
t_temperature = linspace(0, ...
    milliseconds(tempData_img.time(end)-tempData_img.time(1)), ...
    height(tempData_img));

% change time axes for voltage and image to standard vector, do not use
% time, it is not possible to align all of them together (different precisions)

% get stimulus during image recording
m_stimulus_img = ...
    stimulus_img_t<(imgData.pxSzT*size(imgData.wholeImgFluoXT,2)+milliseconds(trigg_delay)) & ...
    stimulus_img_t>milliseconds(trigg_delay);
t_stimAndCaff = stimulus_img_t(m_stimulus_img);
t_stimAndCaff = t_stimAndCaff - t_stimAndCaff(1);
stimAndCaff = table(t_stimAndCaff(:), ...
    stimulusData_V.y(m_stimulus_img), ...
    stimulusData_caff.y(m_stimulus_img), ...
    'VariableNames',{'t', 'stimulus', 'caffeinePerfusion'});
stimAndCaff.Properties.VariableUnits = {'ms', '', ''};

% save temperature, stimulus and perfusion data
imgData.tempData_img = tempData_img;
imgData.stimAndCaff = stimAndCaff;
setappdata(mainFig, 'imgData', imgData)

if showLoadedData
    % show figure with temperature recording and line scan image
    scSz = get(0,'ScreenSize');
    figure('Position',[5 1*scSz(4)/4 scSz(3)-10 3*scSz(4)/4])
    tiledlayout(4,1,"TileSpacing","tight")
    % time limits
    t_lims = [0, imgData.pxSzT*size(imgData.wholeImgFluoXT,2)];
    % plot stimulation protocol
    ax1 = nexttile;
    plot(t_stimAndCaff, ...
        stimulusData_V.y(m_stimulus_img), 'LineWidth',3)
    set(gca, 'FontSize',22, ...
        'XLim', t_lims, 'XTick',[])
    ylabel({"field stimulation","protocol"}, "FontSize",24)
    % plot caffeine local perfusion
    ax2 = nexttile;
    plot(t_stimAndCaff, ...
        stimulusData_caff.y(m_stimulus_img), 'LineWidth',3)
    set(gca, 'FontSize',22, ...
        'XLim', t_lims, 'XTick',[])
    ylabel("caffeine perfusion", "FontSize",24)
    % plot temperature
    ax3 = nexttile;
    plot(t_temperature, tempData_img.temperature, 'LineWidth',3)
    set(gca, 'FontSize',22, ...
        'XLim', t_lims, ...
        'YLim', [floor(min(tempData_img.temperature(:))) ...
        ceil(max(tempData_img.temperature(:)))], ...
        'XTick',[])
    ylabel("temperature (\circC)", "FontSize",24)
    % show image
    ax4 = nexttile;
    imagesc(imgData.wholeImgFluoXT, ...
        'YData',[1 size(imgData.wholeImgFluoXT,1)], ...
        'XData',[0 imgData.pxSzT*size(imgData.wholeImgFluoXT,2)])
    set(gca, 'FontSize',22, 'YTick',[], ...
        'XLim', t_lims)
    xlabel("time (ms)")
    % link axis
    linkprop([ax1.XAxis, ax2.XAxis, ax3.XAxis, ax4.XAxis],'Limits');
end

set(mainFig, "Pointer","arrow")

end