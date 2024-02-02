function saveSparkDetAnalysis(~, ~, mainFig)
% save sparks analysis
% get data and get path where to save results
imgData = getappdata(mainFig,'imgData');
sparkDetection = getappdata(mainFig,'sparkDetection');
analysisType = getappdata(mainFig,'analysisType');
hObjs = getappdata(mainFig,'hObjs'); 

% check if there is something to save
if isempty(sparkDetection)
    return
else
    if ~isfield(sparkDetection,'eventParams')
        return
    end   
end

% get path to file where to save
[~, name, ~] = fileparts(char(imgData.fileName)); 

if isfield(getappdata(mainFig),'batchProcessing')
    if getappdata(mainFig,'batchProcessing')
        file = char([name,'_SparkDetResults','.xls']);
        path = char(imgData.filePath);
    else
        filter = fullfile(char(imgData.filePath), ...
            char([name,'_SparkDetResults','.xls']));      
        [file,path] = uiputfile(filter);
    end
else
    filter = fullfile(char(imgData.filePath), ...
        char([name,'_SparkDetResults','.xls']));      
    [file,path] = uiputfile(filter);
end

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
imgSzT = imgSzPx(2)*pxSzT;        % ms

% create sparks data table
dataSparks = struct2table(sparkDetection.eventParams);
dataSparks.maskOfAcceptedSparks = sparkDetection.maskOfAcceptedSparks(:);
dataSparks.Properties.VariableUnits = ...
    {'deltaF/F0' 'ms' 'ms' 'ms' 'um' ...
     'deltaF/F0*um^3' 'ms' 'ms' 'deltaF/F0*ms' ...
     '' '' '' '' '' ''};

resultSparksParamsNames = cellfun(@(x,y) [x,' (',y,')'], ...
    dataSparks.Properties.VariableNames, ...
    dataSparks.Properties.VariableUnits, ...
    'UniformOutput',0);

resultSparksParams = [resultSparksParamsNames; ...
                      table2cell(dataSparks)];

% save noise of signal in profiles
noiseOfImage = [
    {'noiseEstimate (STD of image without sparks areas, F/F0)'}, ...
    imgData.stdNoise ];

% save position of crop ROI
if ~isfield(imgData,'cropROIpos')
    posOfCropROI = [ [{'position of cropROI'},{''}]; ...
        [{'x(ms)'}; {'y(px)'}; {'w(ms)'}; {'h(px)'}],...
        num2cell([0;0;0;0]) ];
else
    posOfCropROI = [ [{'position of cropROI'},{''}]; ...
        [{'x(ms)'}; {'y(px)'}; {'w(ms)'}; {'h(px)'}],...
        num2cell(imgData.cropROIpos(:)) ];  
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
    [ [{'sparkSelectionCriteria:'},{''} ];
    [{'spAmpl > # (F/F0)'},str2double(hObjs.h_edit_noise.String) ];
    [{'spFDHM > # (ms)'},str2double(hObjs.h_edit_spFDHM.String) ];
    [{'spFWHM > # (um)'},str2double(hObjs.h_edit_spFWHM.String) ];
    [{'smoothingSpan (ms)'},str2double(hObjs.h_smooth_edit.String) ];
    [{'baseline detection sensitivity'},str2double(hObjs.h_bsDet_edit.String) ] ];
    
sparkFreq = {'spark frequency (#sp*100um-1*s-1)',sparkDetection.sparkFreq};
correctedSparkFreq = {'corrected spark frequency (#sp*100um-1*s-1)',sparkDetection.correctedSparkFreq};


%% save data to xls file
% delete xls file if exist
if exist(path_xls, 'file')
    delete(path_xls)
end

% create cell array to write as xls file
cellArr_info1 = [imgPath; imgName; cell(1,2); ...
                 imgSize_px; imgSize; pxSz_x; ...
                 pxSz_t; blank; noiseOfImage; cell(1,2); ...
                 animal; [exp_notes, cell(numel(exp_notes),1)]];

cellArr_info2 = [sparkFreq; correctedSparkFreq; cell(2,2); ...
                 [imgFiltersUsed, cell(numel(imgFiltersUsed),1)]; cell(2,2); ...
                 posOfCropROI; cell(2,2); ...
                 parametersOfSparkDetection; cell(2,2); ...
                 sparkSelectionCriteria];
% get maximum height of cell arrays, pad with empty cells
maxArrHeight = max([size(cellArr_info1,1), ...
    size(cellArr_info2,1), size(resultSparksParams,1)]);
cellArr_info1 = [cellArr_info1; ...
    cell(uint32(maxArrHeight-size(cellArr_info1,1)),size(cellArr_info1,2))];
cellArr_info2 = [cellArr_info2; ...
    cell(uint32(maxArrHeight-size(cellArr_info2,1)),size(cellArr_info2,2))];
resultSparksParams = [resultSparksParams; ...
    cell(uint32(maxArrHeight-size(resultSparksParams,1)),size(resultSparksParams,2))];

cellArrToWrite = [cellArr_info1, cell(maxArrHeight,1), ...
                  cellArr_info2, cell(maxArrHeight,1), ...
                  resultSparksParams];
% save data
writecell(cellArrToWrite, path_xls, 'Sheet', 'info');

%% create final whole image figures 
scRes = get(0,'ScreenSize');

%filter data, whole image
whImg = imgData.wholeImgFluoXT;

% cropped image, filtered and normalized
crImg = imgData.imgDataXTfluoFN;
if isfield(imgData,'cropROIpos')
    cropROIpos = imgData.cropROIpos;
    
    if isfield(imgData,'crop_s_t')
        crop_s_t = round(imgData.crop_s_t/pxSzT); % in px
    else
        crop_s_t = [];
    end
    
else
    cropROIpos = [];
    crop_s_t = [];
end

t_whImg = linspace(0,(size(whImg,2)-1)*pxSzT,size(whImg,2));


%% plot output whole figure
whImg_fig_final = figure('Tag', 'whImg_fig_final');
set(whImg_fig_final, ...
    'Position',[1 1 scRes(3) scRes(4)], ...
    'Visible','on')

whImg_ax = axes('Parent',whImg_fig_final, ...
    'Units','normalized', 'Position',[0.03 0.7 0.94 0.25]);
image(whImg, 'CDataMapping','scaled', 'Parent',whImg_ax, ...
    'XData',[0 size(whImg,2)*pxSzT]);
if isempty(cropROIpos)
    cropROIpos = [whImg_ax.XLim(1) 
                  whImg_ax.YLim(1) 
                  whImg_ax.XLim(2) 
                  whImg_ax.YLim(2)];
end
rectangle('Position',cropROIpos, 'Parent',whImg_ax,...
    'EdgeColor','w', 'LineWidth',5)
colormap(parula(256))
set(whImg_ax,'FontSize',16)
title(whImg_ax,'whole image, filtered')
xlabel(whImg_ax,'t (ms)')
ylabel(whImg_ax,'x (pixels)')

% analyzed image with detected sparks
cImg_ax = axes('Parent',whImg_fig_final,...
    'Units','normalized', 'Position',[0.03 0.35 0.94 0.25]);
yyaxis(cImg_ax,'left')
cImg_ax.YAxis(1).Color = [0 0 0];

image(crImg,'CDataMapping','scaled', 'Parent',cImg_ax,...
    'YData',[0 size(crImg,1)], 'XData',[0 size(crImg,2)*pxSzT]);
ylabel(cImg_ax, 'x (pixels)')
colormap(cImg_ax, parula(256))

clim(cImg_ax, [floor(prctile(crImg(:),1)*10)/10 ...
               ceil(prctile(crImg(:),99.9)) ])
h_cb = colorbar(cImg_ax,'west');
h_cb.Position = [0.95 0.15 0.01 0.15];
h_cb.Label.String = '\DeltaF/F_0';
%h_cb.Label.Position = [0.9 2.5 0];

yyaxis(cImg_ax,'right')
cImg_ax.YAxis(2).Color = [0 0 0];
hImg = image(crImg,'CDataMapping','scaled', 'Parent',cImg_ax,...
    'YData',[0 size(crImg,1)*pxSzX], ...
    'XData',[0 size(crImg,2)*pxSzT]);
delete(hImg)
ylabel(cImg_ax,'x (\mum)')

xlabel(cImg_ax,'t (ms)')
set(cImg_ax,'FontSize',16) 
title(cImg_ax,'cropped image (white rectangle), filtered & normalized')

yyaxis(cImg_ax,'left')
% copyobj(sparkDetection.detectedEventsRec,cImg_ax)
% plot bounding rectangles of sparks
for i = 1:numel(sparkDetection.detectedEventsRec)
    rectangle('Position',sparkDetection.detectedEventsRec(i).Position,...
        'Parent',cImg_ax,...
        'EdgeColor',sparkDetection.detectedEventsRec(i).EdgeColor,...
        'LineWidth',1);
    text(cImg_ax,sparkDetection.detectedEventsRec(i).Position(1),...
        sparkDetection.detectedEventsRec(i).Position(2),num2str(i),...
        'FontSize',8, 'VerticalAlignment','bottom', 'FontWeight','bold',...
        'Color',sparkDetection.detectedEventsRec(i).EdgeColor)
end

%invisible axes for text 
text_ax = axes('Parent',whImg_fig_final,'Units','normalized','Position',[0.03 0.1 0.94 0.2]);
set(text_ax,'Visible','off')

% create text input
notes = hObjs.h_table_notes.Data;
lineNotes = strjoin(cellfun(@(x) [x,' // '], notes,'UniformOutput',0));

ImgDataPath = fullfile(imgData.filePath,imgData.fileName);
ImgDataPath = strrep(ImgDataPath, '_', '\_');

% add average spark params
dataAvrg = num2cell(mean(dataSparks{dataSparks.maskOfAcceptedSparks,(1:6)},1));
if ~isempty(dataAvrg)
    avrgValuesOfSparks = strjoin( cellfun( @(x,y) sprintf('%s = %0.2f // ',x,y), ...
        resultSparksParamsNames(1:6), dataAvrg,'UniformOutput',0 ) );
    % add greek letters
    avrgValuesOfSparks = strrep(avrgValuesOfSparks, 'deltaF/F0', '\Delta/F_0');
    avrgValuesOfSparks = strrep(avrgValuesOfSparks, 'um', '\mum');
else
    avrgValuesOfSparks = 'no sparks detected with selected criteria';
end

avrgValuesOfSparks = [{'average spark parameters (from accepted sparks):'}; ...
    {avrgValuesOfSparks}; ...
    {[sprintf('spark frequency (accepted sparks) = %0.2f',correctedSparkFreq{2}), ' (#sp*100\mum^{-1}*s^{-1})']}];

txtStr =[ {sprintf('animal: %s',hObjs.popUpMenuAnimal.String{hObjs.popUpMenuAnimal.Value}); ...
    sprintf('%s',lineNotes); ...
    sprintf('ImgDataPath: %s',ImgDataPath)}; ...
    avrgValuesOfSparks ];

text(min(get(text_ax,'XLim')),max(get(text_ax,'YLim')),txtStr,...
    'Parent',text_ax,'FontUnits','points','FontSize',16,'VerticalAlignment','top')

% figure of xy image, position of scanning line and TPP point, scale 
if ~isempty(imgData.imgDataXY)
    
    scanLinePos = imgData.scanLinePos;
    xyCellImg = imgData.imgDataXY;
    
    rp = min(scRes(3:4))/max(scRes(3:4));    
    dres = min(scRes(3:4))*2 - max(scRes(3:4));   
    szAx = min(scRes(3:4))-floor(dres/2);
          
    xyCellImg_h = figure('Tag','xyCellImg');
    set(xyCellImg_h,'Position',[1 1 2*szAx szAx])
    ch1_h =  axes('Parent',xyCellImg_h,'Position',[0.05 0.2 0.6*rp 0.6]);
    image(xyCellImg{1,1},'Parent',ch1_h,'CDataMapping','scaled')
    colormap(ch1_h,'jet')
    set(ch1_h,'XTick',[],'YTick',[])
    if ~isempty(scanLinePos)
        line([scanLinePos(1) scanLinePos(2)],[scanLinePos(3) scanLinePos(4)],'Parent',ch1_h,'LineWidth',2,'Color','y',...
            'LineStyle','-');        
    end
   
    % scale 10 um
    line([ch1_h.XLim(2)-5-10/pxSzX ch1_h.XLim(2)-5],[ch1_h.YLim(2)-5 ch1_h.YLim(2)-5],'Parent',ch1_h,'LineWidth',2,'Color','g',...
            'LineStyle','-'); 
    text(ch1_h.XLim(2)-5-5/pxSzX,ch1_h.YLim(2)-15,'10 \mum','HorizontalAlignment','center',...
        'Color','g','FontSize',14,'Parent',ch1_h,'FontWeight','bold')
    try   
    ch2_h =  axes('Parent',xyCellImg_h,'Units','normalized','Position',[0.6/(scRes(3)/scRes(4)) + 0.1 0.2 0.6/(scRes(3)/scRes(4)) 0.6]);
    image((xyCellImg{2,1}),'Parent',ch2_h,'CDataMapping','scaled')
    colormap(ch2_h,'gray')
    set(ch2_h,'XTick',[],'YTick',[])
    if ~isempty(scanLinePos)
        line([scanLinePos(1) scanLinePos(2)],[scanLinePos(3) scanLinePos(4)],'Parent',ch2_h,'LineWidth',2,'Color','y',...
            'LineStyle','-');        
    end
    % scale 10 um
    line([ch1_h.XLim(2)-5-10/pxSzX ch1_h.XLim(2)-5],[ch1_h.YLim(2)-5 ch1_h.YLim(2)-5],'Parent',ch2_h,'LineWidth',2,'Color','g',...
            'LineStyle','-'); 
    text(ch1_h.XLim(2)-5-5/pxSzX,ch1_h.YLim(2)-15,'10 \mum','HorizontalAlignment','center',...
        'Color','g','FontSize',14,'Parent',ch2_h,'FontWeight','bold')
    catch
    end
else
    xyCellImg_h = [];
end


%% save figures

% get handle of imgage parts analysis
h_imgPartsAnalysis = findobj('-regexp','Tag','parts of image');

% handles of figures of ca events, take only accepted sparks 
h_figsCaEvents = findobj('Type','figure','Tag','CaEventParam');
if ~isempty(h_figsCaEvents)
    eventNum = arrayfun(@(x) sscanf(x.Name,'spark from ROI #: %d'), ...
        h_figsCaEvents, 'UniformOutput',1);
    [~,indx_e] = sort(eventNum, 1, 'ascend');
    h_figsCaEvents = h_figsCaEvents(indx_e);
    % take only accepted sparks
    acceptedSparks = dataSparks(dataSparks.maskOfAcceptedSparks,:);
    % close not accepted ones
    h_figsCaEvents_notAccepted = h_figsCaEvents( ...
        ~ismember(eventNum(indx_e), acceptedSparks.sparkROINum));
    h_figsCaEvents = h_figsCaEvents( ...
        ismember(eventNum(indx_e), acceptedSparks.sparkROINum));
    close(h_figsCaEvents_notAccepted)
end

% get handles of all figures to save
if hObjs.check_saveEventsFigs.Value
    if isempty(h_imgPartsAnalysis)
        allHtoSave = [whImg_fig_final; xyCellImg_h; h_figsCaEvents];
    else
        allHtoSave = [whImg_fig_final; h_imgPartsAnalysis; ...
            xyCellImg_h; h_figsCaEvents];
    end
        
else
    if isempty(h_imgPartsAnalysis)
        allHtoSave = [whImg_fig_final; xyCellImg_h];
    else
        allHtoSave = [whImg_fig_final; h_imgPartsAnalysis; xyCellImg_h];
    end
    close(h_figsCaEvents)
end
 
% get path to save figures 
[pathFigs, nameFigs,~] = fileparts(path_xls);
% output file path with selected format
outputFiguresFilePath = sprintf('%s/%s.%s', pathFigs, nameFigs, ...
    hObjs.outputFileFormatRBgroup.SelectedObject.String);
% delete files if exist
if exist(outputFiguresFilePath, 'file')
    delete(sprintf('%s/%s.pdf', pathFigs, nameFigs))
end
% ask if delete also other file formats if present
if exist(sprintf('%s/%s.gif', pathFigs, nameFigs), 'file') || ...
        exist(sprintf('%s/%s.tif', pathFigs, nameFigs), 'file') || ...
        exist(sprintf('%s/%s.ps', pathFigs, nameFigs), 'file')
    answer = questdlg( ...
        'Would you like to delete analysis file with different format (.pdf, .gif, .ps)?', ...
	    'Output analysis files in different file format.', ...
	    'YES', 'NO', 'NO');
    switch answer
        case 'YES'
            delete(sprintf('%s/%s.gif', pathFigs, nameFigs))
        case 'NO'
    end
end

% save figures as .gif or .pdf
if ~isempty(allHtoSave)
    for i = 1:numel(allHtoSave)     
        % set upfigures for printing
        set(allHtoSave(i), 'Units','pixels');
        posFig = get(allHtoSave(i),'Position');
        set(allHtoSave(i),'PaperPositionMode','Auto',...
                          'PaperUnits','points',...
                          'PaperSize',[posFig(3), posFig(4)])   
        % create file 
        if ~exist(outputFiguresFilePath, 'file')
            exportgraphics(allHtoSave(i), ...
                outputFiguresFilePath, ...
                'Resolution',str2double(hObjs.h_edit_res.String))
        else
            % append
            exportgraphics(allHtoSave(i), ...
                outputFiguresFilePath, ...
                'Append',true, ...
                'Resolution',str2double(hObjs.h_edit_res.String))
        end                                    
    end  
    close(allHtoSave)
end

%%
close(findobj('Type','figure','Name','Experiment prewiev'));

% clear window
if isfield(hObjs, 'ax_img_sparks')
    cla(hObjs.ax_img_sparks);
end
cla(hObjs.ax_img);
cla(hObjs.h_ax_transCh);

delete(get(hObjs.ax_img,'Children'))
cla(hObjs.ax_prof);

set(mainFig,'Pointer','arrow')
drawnow

% % if this is re-analyzed image, called from spark recovery analysis than close
% % window
% if isfield(getappdata(mainFig),'reAnalyzeImg')
%     if getappdata(mainFig,'reAnalyzeImg') == true
%         close(mainFig)
%     end
% end

end






