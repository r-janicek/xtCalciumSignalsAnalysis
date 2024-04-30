function saveAnalysis(~,~,mainFig)
% save data to xls file, and crete pdf with results

%% get data and path of file where to save
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

selectedROIs = getappdata(mainFig,'selectedROIs');

if ~isfield(getappdata(mainFig),'selectedROIs')
    return
else
    % check if there is something to save
    if isempty(selectedROIs)
        return
        
    else
        if any( cellfun(@(x) isempty(x), selectedROIs.analysis ))
            warndlg('Analyze all selected events!')
            return
        end
    end
end

% get path to file where to save
[~,name,~] = fileparts(char(imgData.fileName)); 

filter = fullfile(char(imgData.filePath), ...
    char([name,'_GlobEvntsResults','.xls']));      
[file,path] = uiputfile(filter);

if path == 0
    return
end 

set(mainFig,'Pointer','watch')
drawnow

% set up path for xls file
path_xls = fullfile(path, file);

%% get all data and prepare them to save in xls file 
imgSzPx = size(imgData.imgDataXTfluoFN);

filePath = imgData.filePath;
fileName = imgData.fileName;

pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;

imgSzX = imgSzPx(1)*pxSzX;        % um
imgSzT = imgSzPx(2)*pxSzT-pxSzT;  % ms

t = imgData.t;
t = t(:);

% create table with all results
for i = 1:height(selectedROIs)   
  
    % get position of ROI
    posOfROIs = selectedROIs.positionOfRoi{i}; % y in pixels; x in ms
    nameOfROI = selectedROIs.roiName{i};
    
    % get analysis of events
    eventAnalysis = selectedROIs.analysis{i};
    type = eventAnalysis.type;
    t = selectedROIs.dataROIs(i).t;
    switch type
        case 'wave'
            waveSpeed = eventAnalysis.waveAnalysis.waveSpeed;
            prof = eventAnalysis.waveAnalysis.deskewedWaveProf;
            splFit = eventAnalysis.waveAnalysis.fitOfEvent.splFit;
            calcParamsFromFit = eventAnalysis.waveAnalysis.fitOfEvent.calcParamsFromFit;
            
        case 'caffeine'
            waveSpeed = [];
            prof = eventAnalysis.caffeineAnalysis.deskewedCaffeineProf;
            splFit = eventAnalysis.caffeineAnalysis.fitOfEvent.splFit;
            calcParamsFromFit = eventAnalysis.caffeineAnalysis.fitOfEvent.calcParamsFromFit;
               
        case 'transient'
            if isempty(regexp(selectedROIs.roiName{i},'spon','once'))
                type = 'transient_trigg';
            else
                type = 'transient_spon';
            end
            waveSpeed = [];
            prof = eventAnalysis.transientAnalysis.transientsProf;
            splFit = eventAnalysis.transientAnalysis.fitOfEvent.splFit;
            calcParamsFromFit = eventAnalysis.transientAnalysis.fitOfEvent.calcParamsFromFit;      
    end
    
    % estimation of noise in F/F0 units
    noiseOfProf = {std( prof( prof<(mean(prof)+1*std(prof)) ) )};  
    
    % save profile of event
    profiles{i}(:,1) = [{sprintf('%s--t (ms)',nameOfROI)}; ...
        num2cell(t(:))];
    profiles{i}(:,2) = [{sprintf('%s--tProfile (deltaF/F0)',nameOfROI)}; ...
        num2cell(prof(:))];
    profiles{i}(:,3) = [{sprintf('%s--tProfileFit (deltaF/F0)',nameOfROI)}; ...
        num2cell(splFit(:))];
    
    % put in cell array
    type = {type};
    nameOfROI = {nameOfROI};
    waveSpeed = {waveSpeed};
    
    % if there is more analyzed peaks per event 
    if size(calcParamsFromFit,1) > 2
        type = repmat(type,size(calcParamsFromFit,1)-1,1);
        waveSpeed = repmat(waveSpeed,size(calcParamsFromFit,1)-1,1);
        posOfROIs = repmat(posOfROIs,size(calcParamsFromFit,1)-1,1);
        nameOfROI = repmat(nameOfROI,size(calcParamsFromFit,1)-1,1);
        noiseOfProf = repmat(noiseOfProf,size(calcParamsFromFit,1)-1,1);
    end

    % recalculate time positions in calcParamsFromFit, so they are starting
    % from the beginning of image (cropped)
    calcParamsFromFit(2:end,[1,2,4,5,6,17]) = ...
        num2cell(cell2mat(calcParamsFromFit(2:end,[1,2,4,5,6,17])) + posOfROIs(1));
    
    % create results array
    if exist('result','var')
        
        result = [ result; ...
            [ num2cell(posOfROIs), nameOfROI, type, noiseOfProf, waveSpeed,...
           calcParamsFromFit(2:end,1:end-1) ] ];   
    else
        nn = {'posOfROI_x (ms)', 'posOfROI_y (px)', 'posOfROI_w (ms)', 'posOfROI_h (px)', ...
            'nameOfROI', 'type', 'noiseEstimate (STD of profile, deltaF/F0)', ...
            'waveSpeed (um/s)'};
        nnCalcParams = calcParamsFromFit(1,1:end-1);
        
        result = [ [nn,nnCalcParams]; ...
           [ num2cell(posOfROIs), nameOfROI, type, noiseOfProf, waveSpeed,...
           calcParamsFromFit(2:end,1:end-1) ] ];  
    end
end

% save position of crop ROI
if isfield(imgData,'cropROIpos')
    posOfCropROI = [{'position of cropROI'};{'x(ms); y(px); w(ms); h(px)'};...
        num2cell(imgData.cropROIpos(:))];
else
    posOfCropROI = [{'position of cropROI'};{'x(ms); y(px); w(ms); h(px)'};...
        num2cell(nan(4,1))];
end

% data animal and notes
animal = hObjs.popUpMenuAnimal.String(hObjs.popUpMenuAnimal.Value);

exp_notes = hObjs.h_table_notes.Data;
exp_notes = [{'exp_notes:'};exp_notes];
animal = [{'animal:'},animal];

% filters used on image
imgFiltersUsed = imgData.imgFiltersUsed;
imgFiltersUsed = [{'filters applied to image:'};{imgFiltersUsed}];

imgPath = {'imgPath:',filePath};
imgName = [{'imgName:'},fileName];
pxSz_x = {'px size x (um)',pxSzX};
pxSz_t = {'px size t (ms)',pxSzT};

imgSize = [{'image size x (um)',imgSzX};{'image size t (ms)',imgSzT}];
imgSize_px = [{'image size x (pixels)',imgSzPx(1)};{'image size t (pixels)',imgSzPx(2)}];

blank = {'blank:',imgData.blank};


%% save data to xls file
% delete xls file if exist
if exist(path_xls, 'file')
    delete(path_xls)
end
% create cell array (INFO) to write as xls file
cellArr_info1 = [imgPath; imgName; cell(1,2); ...
                 imgSize_px; imgSize; pxSz_x; ...
                 pxSz_t; blank; ...
                 animal; [exp_notes, cell(numel(exp_notes),1)]];

cellArr_info2 = [[imgFiltersUsed, cell(numel(imgFiltersUsed),1)]; ...
                 cell(2,2); ...
                 [posOfCropROI, cell(numel(posOfCropROI),1)]; ...
                 cell(2,2)];
% get maximum height of cell arrays, pad with empty cells
maxArrHeight = max([size(cellArr_info1,1), ...
    size(cellArr_info2,1), size(result,1)]);
cellArr_info1 = [cellArr_info1; ...
    cell(uint32(maxArrHeight-size(cellArr_info1,1)),size(cellArr_info1,2))];
cellArr_info2 = [cellArr_info2; ...
    cell(uint32(maxArrHeight-size(cellArr_info2,1)),size(cellArr_info2,2))];
result = [result; ...
    cell(uint32(maxArrHeight-size(result,1)),size(result,2))];

cellArrToWrite_info = [cellArr_info1, cell(maxArrHeight,1), ...
                  cellArr_info2, cell(maxArrHeight,1), ...
                  result];
% save info data
writecell(cellArrToWrite_info, path_xls, 'Sheet','info');

% check whether save profiles and fits
if hObjs.check_saveProfsAndFits.Value
    % prepare cell array with all profiles
    maxSz = max(cellfun(@(x) size(x,1),profiles));
    allProfilesData = [];
    for i = 1:numel(profiles)
        allProfilesData = [allProfilesData, ... 
            [ [profiles{i}; ...
                cell(uint32(maxSz-size(profiles{i},1)),size(profiles{i},2))], ...
                cell(uint32(maxSz),1) ] ];
    end
    % save all profiles data
    writecell(allProfilesData, path_xls, 'Sheet','profilesOfEvents');    
end

% save temperature recorded during experiment
if hObjs.check_addTemperature.Value
    % check if exists
    if isfield(imgData, 'tempData_img')
        writetable(imgData.tempData_img, path_xls, 'Sheet','temperature');
    end
end


%% create final whole image figures 
scRes = get(0,'ScreenSize');

%filter data, whole image
whImg = imgData.wholeImgFluoXT;

% cropped image, filtered and normalized
crImg = imgData.imgDataXTfluoFN;
if isfield(imgData,'cropROIpos')
    cropROIpos = imgData.cropROIpos;
    
    if isfield(imgData,'crop_s_t')
        crop_s_t = imgData.crop_s_t; %round(imgData.crop_s_t/pxSzT); % in px
    else
        crop_s_t = 0;
    end    
else
    cropROIpos = [];
    crop_s_t = 0;
end
t_whImg = linspace(0,(size(whImg,2)-1)*pxSzT,size(whImg,2));
 

%% plot output whole figure
whImg_fig_final = figure('Tag','whImg_fig_final');
set(whImg_fig_final,'Position',[1 1 scRes(3) scRes(4)])

% plot whole image with white rectangle where it was cropped

% add temperature recorded during experiment
if hObjs.check_addTemperature.Value
    % check if temperature recording exists
    if isfield(imgData, 'tempData_img')
        temp_ax = axes('Parent',whImg_fig_final, ...
            'Units','normalized', ...
            'Position',[0.04 0.895 0.94 0.09]);
        plot(imgData.tempData_img.time, ...
            imgData.tempData_img.temperature, ...
            'Parent',temp_ax, 'Color','b', ...
            'LineWidth',2);
        ylabel(temp_ax, 'temperature (\circC)', 'FontSize',14)
        set(temp_ax, 'FontSize',14)
        set(temp_ax, 'XTick',[], ...
            'XLim',[min(imgData.tempData_img.time) ...
            max(imgData.tempData_img.time)],...
            'YLim',[floor(min(imgData.tempData_img.temperature)) ...
            ceil(max(imgData.tempData_img.temperature))])
    end
end

whImg_ax = axes('Parent',whImg_fig_final, 'Units','normalized', ...
    'Position',[0.04 0.66 0.94 0.22]);
image(whImg, 'CDataMapping','scaled', 'Parent',whImg_ax, ...
    'XData',[0 size(whImg,2)*pxSzT]);
if isempty(cropROIpos)
    cropROIpos = [whImg_ax.XLim(1), whImg_ax.YLim(1), ...
        whImg_ax.XLim(2), whImg_ax.YLim(2)];
end
rectangle('Position',cropROIpos, 'Parent',whImg_ax,...
    'EdgeColor','w', 'LineWidth',5)
colormap(parula(256))
set(whImg_ax, 'XTick',[], 'FontSize',14)
title('whole image, filtered')
%xlabel('t (ms)')
ylabel('x (pixels)')

% add stimulus and caffeine perfusion
if hObjs.check_addTemperature.Value
    if isfield(imgData, 'stimAndCaff')
        stim_ax = axes('Parent',whImg_fig_final, ...
            'Units','normalized', ...
            'Position',[0.04 0.58 0.94 0.06]);
        plot(imgData.stimAndCaff.t-crop_s_t, ...
            imgData.stimAndCaff.stimulus, ...
            'Parent',stim_ax, 'Color','k', ...
            'LineWidth',2);
        ylabel(stim_ax, 'stimulus', 'FontSize',14)
        set(stim_ax, 'FontSize',14)
        set(stim_ax, 'XTick',[], 'YTick',[],  ...
            'XLim',[imgData.t(1) imgData.t(end)],...
            'YLim',[floor(min(imgData.stimAndCaff.stimulus)) ...
            ceil(max(imgData.stimAndCaff.stimulus))])
        caff_ax = axes('Parent',whImg_fig_final, ...
            'Units','normalized', ...
            'Position',[0.04 0.515 0.94 0.06]);
        plot(imgData.stimAndCaff.t-crop_s_t, ...
            imgData.stimAndCaff.caffeinePerfusion, ...
            'Parent',caff_ax, 'Color','k', ...
            'LineWidth',2);
        ylabel(caff_ax, 'caffeine', 'FontSize',14)
        set(caff_ax, 'FontSize',14)
        set(caff_ax, 'XTick',[], 'YTick',[],...
            'XLim',[imgData.t(1) imgData.t(end)],...
            'YLim',[floor(min(imgData.stimAndCaff.caffeinePerfusion)) ...
            ceil(max(imgData.stimAndCaff.caffeinePerfusion))])
        move_up_axes = 0; 
    else
        move_up_axes = 0.135; 
    end
else
    move_up_axes = 0.135;
end

% make sure it is showing whole image
hObjs.ax_img.XLim = [imgData.t(1) imgData.t(end)];
% copy filtered and normalized cropped image
ax_CroppedImg = copyobj(hObjs.ax_img,whImg_fig_final);
ax_CroppedImg.Position = [0.04 0.275+move_up_axes 0.94 0.22];
title(ax_CroppedImg, ...
    'cropped image (white rectangle), filtered & normalized')
% add colorbar
colormap(ax_CroppedImg, parula(256))
h_CroppedImg = findall(ax_CroppedImg, 'Type','Image');
clim(ax_CroppedImg, getAxisLimits(h_CroppedImg.CData, 1))
h_cb = colorbar(ax_CroppedImg,'south');
h_cb.Position = [0.92 0.035+move_up_axes 0.06 0.01];
h_cb.Label.String = '\DeltaF/F_0';

% copy profile
ax_imgProf = copyobj(hObjs.ax_prof,whImg_fig_final);
ax_imgProf.Position = [0.04 0.10+move_up_axes 0.94 0.17];

% show time to first wave from beginning of last triggered transient
m_triggTrans = strcmp(result(:,6), 'transient_trigg' );
m_wave = strcmp(result(:,6), 'wave' );

transRes = cell2mat(result(2:end,9));
transRes = [nan;transRes];
transRes(~m_triggTrans) = nan;

[lastTrans_t0,indTr]= max(transRes);
firstWave_t0 = min([result{m_wave,9}]);

try
    yL = result{indTr,11} + result{indTr,16};
    firstWaveLat = firstWave_t0 - lastTrans_t0;
    
    line(ax_imgProf,'XData',[lastTrans_t0 lastTrans_t0],...
        'YData',[ax_imgProf.YLim(1) yL],...
        'Color','r','LineWidth',1,'LineStyle',':')
    line(ax_imgProf,'XData',[firstWave_t0 firstWave_t0],...
        'YData',[ax_imgProf.YLim(1) yL],...
        'Color','r','LineWidth',1,'LineStyle',':')
    line(ax_imgProf,'XData',[lastTrans_t0 firstWave_t0],...
        'YData',[yL yL],'Color','r','LineWidth',1)
    text(ax_imgProf,lastTrans_t0+firstWaveLat/2,yL,...
        sprintf('1.wave latency = %g (ms)',firstWaveLat),...
        'VerticalAlignment','top',...
        'HorizontalAlignment','center',...
        'FontSize',14)
catch
    
end

% invisible axes for text 
text_ax = axes('Parent',whImg_fig_final, ...
    'Units','normalized', 'Position',[0.03 0.01+move_up_axes 0.94 0.07]);
set(text_ax,'Visible','off')

notes = hObjs.h_table_notes.Data;
lineNotes = strjoin(cellfun(@(x) [x,' // '], notes,'UniformOutput',0));

ImgDataPath = fullfile(imgData.filePath,imgData.fileName);
ImgDataPath = strrep(ImgDataPath, '_', '\_');

txtStr = {
    sprintf('animal: %s',hObjs.popUpMenuAnimal.String{hObjs.popUpMenuAnimal.Value});...
    sprintf('%s',lineNotes);...
    sprintf('ImgDataPath: %s',ImgDataPath)};

text(min(get(text_ax,'XLim')),max(get(text_ax,'YLim')),txtStr,...
    'Parent',text_ax,'FontUnits','points','FontSize',16,...
    'VerticalAlignment','top')

%% figure of xy image, position of scanning line and TPP point, scale
if ~isempty(imgData.imgDataXY)
    
    scanLinePos = imgData.scanLinePos;
    xyCellImg = imgData.imgDataXY;
    
    rp = min(scRes(3:4))/max(scRes(3:4));    
    dres = min(scRes(3:4))*2 - max(scRes(3:4));   
    szAx = min(scRes(3:4))-floor(dres/2);
          
    xyCellImg_h = figure('Tag','xyCellImg');
    set(xyCellImg_h, 'Position',[1 1 2*szAx szAx])
    ch1_h =  axes('Parent',xyCellImg_h, 'Position',[0.05 0.2 0.6*rp 0.6]);
    image(xyCellImg{1,1}, 'Parent',ch1_h, 'CDataMapping','scaled')
    colormap(ch1_h,'jet')
    set(ch1_h,'XTick',[],'YTick',[])
    if ~isempty(scanLinePos)
        line([scanLinePos(1) scanLinePos(2)], ...
            [scanLinePos(3) scanLinePos(4)], ...
            'Parent',ch1_h, 'LineWidth',2, 'Color','y', 'LineStyle','-');        
    end

    % scale 10 um
    line([ch1_h.XLim(2)-5-10/pxSzX ch1_h.XLim(2)-5], ...
        [ch1_h.YLim(2)-5 ch1_h.YLim(2)-5], ...
        'Parent',ch1_h, 'LineWidth',2, 'Color','g', 'LineStyle','-'); 
    text(ch1_h.XLim(2)-5-5/pxSzX, ch1_h.YLim(2)-15, ...
        '10 \mum', 'HorizontalAlignment','center',...
        'Color','g', 'FontSize',12, 'Parent',ch1_h, 'FontWeight','bold')
    try   
    ch2_h =  axes('Parent',xyCellImg_h, ...
        'Units','normalized', ...
        'Position',[0.6/(scRes(3)/scRes(4))+0.1 0.2 0.6/(scRes(3)/scRes(4)) 0.6]);
    image((xyCellImg{2,1}), 'Parent',ch2_h, 'CDataMapping','scaled')
    colormap(ch2_h,'gray')
    set(ch2_h, 'XTick',[], 'YTick',[])
    if ~isempty(scanLinePos)
        line([scanLinePos(1) scanLinePos(2)], ...
            [scanLinePos(3) scanLinePos(4)], ...
            'Parent',ch2_h, 'LineWidth',2, 'Color','y', 'LineStyle','-');        
    end
    % scale 10 um
    line([ch1_h.XLim(2)-5-10/pxSzX ch1_h.XLim(2)-5], ...
        [ch1_h.YLim(2)-5 ch1_h.YLim(2)-5], ...
        'Parent',ch2_h, 'LineWidth',2, 'Color','g', 'LineStyle','-'); 
    text(ch1_h.XLim(2)-5-5/pxSzX, ch1_h.YLim(2)-15, ...
        '10 \mum', 'HorizontalAlignment','center',...
        'Color','g', 'FontSize',12, 'Parent',ch2_h, 'FontWeight','bold')
    catch
    end
else
    xyCellImg_h = [];
end


%% save figures

% get handles of figures of individial events analysis
h_figsEventsAnalysis = findall(0, 'Type','figure', ...
    'Tag','imgsOfEventsFromAnalysis');
% make them visible again
arrayfun(@(x) set(x,'Visible','on'), h_figsEventsAnalysis)

% get handles of all figures to save
if hObjs.check_saveEventsFigs.Value
    if isempty(h_figsEventsAnalysis)
        allHtoSave = [whImg_fig_final; xyCellImg_h];
    else
        allHtoSave = [whImg_fig_final; h_figsEventsAnalysis; xyCellImg_h];
    end    
else   
    allHtoSave = [whImg_fig_final;xyCellImg_h];
end
           
% get path to save figures 
[pathFigs, nameFigs,~] = fileparts(path_xls);

% delete files if exist
if exist(sprintf('%s/%s.pdf',pathFigs,nameFigs), 'file')
    delete(sprintf('%s/%s.pdf',pathFigs,nameFigs))
end

if exist(sprintf('%s/%s.ps',pathFigs,nameFigs), 'file')
    delete(sprintf('%s/%s.ps',pathFigs,nameFigs))
end

% save figures as .ps and also .pdf
if ~isempty(allHtoSave)
    for i = 1:numel(allHtoSave)      
        % set upfigures for printing
        set(allHtoSave(i),'Units','pixels');
        posFig = get(allHtoSave(i),'Position');
        set(allHtoSave(i),'PaperPositionMode','Auto',...
                          'PaperUnits','points',...
                          'PaperSize',[posFig(3), posFig(4)])
        
        % using print option and postscript     
        % create file 
        if ~exist(sprintf('%s/%s.pdf',pathFigs,nameFigs), 'file')  
            % print(allHtoSave(i), ...
            %      sprintf('%s/%s.ps',pathFigs,nameFigs), '-dpsc')
            exportgraphics(allHtoSave(i), ...
                sprintf('%s/%s.pdf',pathFigs,nameFigs), ...
                'Resolution',str2double(hObjs.h_edit_res.String))
        else
             % print(allHtoSave(i), '-dpsc', '-append', ...
             %     sprintf('%s/%s.ps',pathFigs,nameFigs))
            exportgraphics(allHtoSave(i), ...
                sprintf('%s/%s.pdf',pathFigs,nameFigs), ...
                'Resolution',str2double(hObjs.h_edit_res.String), ...
                'Append',true)
        end                                                                       
    end   

    close(allHtoSave)
end

%%
close(findobj('Type','figure', 'Name','analysis of selected ROIs'));

% clear window
cla(hObjs.ax_img);
cla(hObjs.ax_prof);

delete(get(hObjs.ax_img, 'Children'))
cla(hObjs.ax_prof);

set(mainFig,'Pointer','arrow')
drawnow

% if this re-analyzed image, called from spark recovery analysis than close
% window
if isfield(getappdata(mainFig),'reAnalyzeImg')
    if getappdata(mainFig,'reAnalyzeImg') == true
        close(mainFig)
    else
        % remove analysis from application data
        rmappdata(mainFig,'selectedROIs')
    end
else
    % remove analysis from application data
    rmappdata(mainFig,'selectedROIs')
end

end

