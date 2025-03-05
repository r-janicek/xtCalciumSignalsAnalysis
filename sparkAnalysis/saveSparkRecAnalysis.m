function saveSparkRecAnalysis(~,~,mainFig)
% save data from recovery analysis
imgData = getappdata(mainFig,'imgData');
profileAnalysis = getappdata(mainFig,'profileAnalysis');
analysisType = getappdata(mainFig,'analysisType');
hObjs = getappdata(mainFig,'hObjs'); 

% check if there is something to save
if isempty(profileAnalysis)
    return  
else
    if ~isfield(profileAnalysis,'selectedROIs')
        return 
    else
        selectedROIs = profileAnalysis.selectedROIs;
        if ~(any(strcmp(selectedROIs.Properties.VariableNames, ...
                'AnalysisResult'))) %|| ...
              %any(strcmp(selectedROIs.Properties.VariableNames,'imgPartsAnalysis')))
            return
        end
    end   
end

% get path to file where to save
[~, name,~] = fileparts(char(imgData.fileName)); 
filter = fullfile(char(imgData.filePath), ...
    char([name,'_SparkRecResults','.xls']));      
[file, path] = uiputfile(filter);

if path == 0
    return
end 

set(mainFig,'Pointer','watch')
drawnow
% set up path for xls file
path_xls = fullfile(path,file);

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

for i = 1:height(selectedROIs)   
  
    posOfROIs = repmat(selectedROIs.positionOfRoi(i),height(selectedROIs.AnalysisResult{i}),1);
    posOfROIs_1 = repmat(selectedROIs.positionOfRoi(i),height(selectedROIs.allPeaksData{i}),1);
    
    pp = selectedROIs.normProf{i};
    
    profiles(:,i) = [{sprintf('x_pos %g (px, deltaF/F0)',selectedROIs.positionOfRoi(i))}; ...
                      num2cell(pp)];
    
    noiseOfProf(i,1) = {selectedROIs.positionOfRoi(i)};
    noiseOfProf(i,2) = {std( pp( pp<(mean(pp)+1.5*std(pp)) ) )};  % noise in F/F0 units
                  
    if exist('result','var')
        result = [result; ...
            num2cell([posOfROIs, selectedROIs.AnalysisResult{i}{:,:}])];
        result_allPeaksData = [result_allPeaksData; ...
            num2cell([posOfROIs_1, selectedROIs.allPeaksData{i}{:,:}])];
    else
        nn = {'posOfROI (x,pixels)', 'amplitudeRatio (A2/A1)', ...
              'sparkToSparkDelay (ms)', 'acceptedPair',...
              'sparkAmplRatio (S2/S1)', 'sparkMassRatio (S2/S1)',...
              'amplRatioFit','sparkToSparkDelayFit (ms)', ...
              'TTPratioFit', 'FDHMratioFit', 'AUCratioFit', 'tauDratioFit',...
              'sparkToSparkDelayDerivation', 'amplDerivationRatio', ...
              'sparkToSparkDelayDerivationFit', 'amplDerivationRatioFit'}; 
        
        nn_d = {'posOfROI (x,pixels)', 'baseline (deltaF/F0)', ...
                'amplitude (deltaF/F0)', 'peakPos (ms)', ...
                'sparkAmplitude (deltaF/F0)', 'sparkTTP (ms)', ...
                'sparkTTP_fittedLine (ms)', ...
                'sparkFDHM (ms)', 'sparkFWHM (um)', ...
                'sparkMass (deltaF/F0*um^3)', ...
                'sparkTauD (ms)', 'sparkAUC (um*ms*deltaF/F0)', ...
                'amplitudeFit (deltaF/F0)', 'peakPosFit (ms)', ...
                'TTPfit (ms)', 'FDHMfit (ms)', 'AUC (ms*deltaF/F0)',...
                'tauDfit (ms)',...
                'firstDerivMaxVal (deltaF/F0)', 'firstDerivMaxValPos (ms)', 'firstDerivMaxValFit (deltaF/F0)', 'firstDerivMaxValPosFit (ms)'};
                      
        result = [nn; ...
            num2cell([posOfROIs, selectedROIs.AnalysisResult{i}{:,:}])];
        result_allPeaksData = [nn_d; ...
            num2cell([posOfROIs_1, selectedROIs.allPeaksData{i}{:,:}])];   
    end
    clearvars posOfROIs posOfROIs_1 pp   
end

% save noise of signal in profiles
noiseOfProf = [{'posOfROI (x,pixels)', ...
    'noiseEstimate (STD of profile, F/F0)'}; noiseOfProf];
% save position of crop ROI
posOfCropROI = [{'position of cropROI'};{'x(ms); y(px); w(ms); h(px)'};...
    num2cell(imgData.cropROIpos(:))];
% data animal and notes
animal = hObjs.popUpMenuAnimal.String(hObjs.popUpMenuAnimal.Value);

exp_notes = hObjs.h_table_notes.Data;
exp_notes = [{'exp_notes:'};exp_notes];
animal = [{'animal:'},animal];

% selected profiles data
profiles = [[{'t (ms)'}; num2cell(t)],profiles];

% filters used on image
imgFiltersUsed = imgData.imgFiltersUsed;
imgFiltersUsed = [{'filters applied to image:'};{imgFiltersUsed}];

imgPath = {'imgPath:',filePath};
imgName = [{'imgName:'},fileName];
pxSz_x = {'px size x (um)',pxSzX};
pxSz_t = {'px size t (ms)',pxSzT};

imgSize = [{'image size x (um)',imgSzX}; ...
           {'image size t (ms)',imgSzT}];
imgSize_px = [{'image size x (pixels)',imgSzPx(1)}; ...
              {'image size t (pixels)',imgSzPx(2)}];
blank = {'blank:',imgData.blank};

diameterAroundCentre = [ ...
    {'diameter around spark centre (um):',str2double(hObjs.h_edit_averageWidth.String)};...
    {'# of pixel to calculate average profile:',profileAnalysis.numOfPxAvrg} ];

%% save electrophysiology data, only when data preview is used
if isfield(getappdata(mainFig), 'electroPhys')
    
    electroPhys = getappdata(mainFig,'electroPhys');
    delayImgToCurr = electroPhys.delayImgToCurr; % in ms

    % electrophysiology data 
    currentAmpl = electroPhys.currentAmpl;
        
    if isfield(electroPhys,'cellParameters')
        cellParameters = electroPhys.cellParameters;
        currentDensity = cell2mat(currentAmpl(2:end,3))./cellParameters{1,1}; % pA/pF
        electro = [ [currentAmpl,[{'currDensity (pA/pF)'};num2cell(currentDensity)]];...
                    [{'average'},{nan},num2cell(mean([currentAmpl{2:end,3}])),num2cell(mean(currentDensity))] ];
    else
        electro = [currentAmpl;...
                    [{'average'},{nan},num2cell(mean([currentAmpl{2:end,3}]))]];
    end      
    delayImg = {'delayOfImgToCurrent (ms)',delayImgToCurr};
else
    delayImg = cell(1,2);
end


%% save data to xls file
% delete xls file if exist
if exist(path_xls, 'file')
    delete(path_xls)
end
% create cell array (INFO) to write as xls file
cellArr_info1 = [imgPath; imgName; cell(1,2); ...
                 imgSize_px; imgSize; pxSz_x; ...
                 pxSz_t; blank; diameterAroundCentre; ...
                 delayImg; cell(1,2); ...
                 animal; [exp_notes, cell(numel(exp_notes),1)]];

cellArr_info2 = [[imgFiltersUsed, cell(numel(imgFiltersUsed),1)]; ...
                 cell(2,2); ...
                 [posOfCropROI, cell(numel(posOfCropROI),1)]; ...
                 cell(2,2); ...
                 noiseOfProf; cell(2,2)];
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

% save all detected peaks data
writecell(result_allPeaksData, path_xls, 'Sheet','allPeaksData');

% check whether save profiles and fits
if hObjs.check_saveProfsAndFits.Value
    % save all profiles data
    writecell(profiles, path_xls, 'Sheet','profiles');
    % save all fitted events data
    try
        if any(strcmp(selectedROIs.Properties.VariableNames, ...
                'wholeProfileFit'))
            for i = 1:height(selectedROIs)
                fit = selectedROIs.wholeProfileFit{i}.profFit;
                
                fit_data = [[{'time (ms)'}; num2cell(fit.t)],...
                    [{'wholeFit'}; num2cell(fit.wholeFit)],...
                    [{'allEventsFit'}; num2cell(fit.allEventsFit)],...
                    [{'baselineFit'}; num2cell(fit.baselineFit)],...
                    fit.individualEventsFits];
                
                sheetName = sprintf('fitDataROIpos(%d px)', ...
                    selectedROIs.positionOfRoi(i));
                
                writecell(fit_data, path_xls, 'Sheet',sheetName);
                
                clearvars fit_data
            end
        end    
    catch
    end  
end

% save current ...
if exist('electro','var')
    writecell(electro, path_xls, 'Sheet','electroPhys');
end

if exist('cellParameters','var')
    writecell([cellParameters.Properties.VariableNames;...
        cellParameters.Properties.VariableUnits;...
        table2cell(cellParameters)], path_xls, ...
        'Sheet','currentAndCell');  
end
clearvars delayImgToCurr xt_cur V I electrophys delayImg


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
whImg_fig_final = figure('Tag','whImg_fig_final');
set(whImg_fig_final,'Position',[1 1 scRes(3) scRes(4)])
whImg_axFontSz = 10;

if isfield(getappdata(mainFig),'electroPhys') 
    electroPhys = getappdata(mainFig,'electroPhys');
    IcaLimits = electroPhys.IcaLimits;   
    delayImgToCurr = electroPhys.delayImgToCurr; % in ms
    xt_cur = electroPhys.xt_cur;
    V = electroPhys.voltage;
    ICa = electroPhys.current;
    [~, pos_delay] = min(abs(xt_cur-delayImgToCurr));
       
    vol_ax = axes('Parent',whImg_fig_final, ...
        'Units','normalized', ...
        'Position',[0.03 0.93 0.94 0.06]);
    plot(xt_cur(pos_delay:end), V(pos_delay:end), ...
        'Parent',vol_ax, 'Color','k', 'LineWidth',2);
    ylabel(vol_ax, 'voltage (mV)')
    set(vol_ax, 'XTick',[], 'FontSize',whImg_axFontSz)
    set(vol_ax, 'XLim',[xt_cur(pos_delay) ...
        max(size(whImg,2)*pxSzT+delayImgToCurr,xt_cur(end))],...
        'YLim',[min(V(pos_delay:end))-5 max(V(pos_delay:end))+5])
    
    curr_ax = axes('Parent',whImg_fig_final, ...
        'Units','normalized', ...
        'Position',[0.03 0.795 0.94 0.13]);
    plot(xt_cur(pos_delay:end), ICa(pos_delay:end), ...
        'Parent',curr_ax, 'Color','r');
    ylabel(curr_ax, 'current (pA)')
    set(curr_ax, 'XTick',[], 'FontSize',whImg_axFontSz)
    set(curr_ax, 'XLim',[xt_cur(pos_delay) ...
        max(size(whImg,2)*pxSzT+delayImgToCurr,xt_cur(end))],...
        'YLim',IcaLimits)
    
    whImg_ax = axes('Parent',whImg_fig_final, ...
        'Units','normalized', ...
        'Position',[0.03 0.55 0.94 0.22]);
    image(whImg, 'CDataMapping','scaled', ...
        'Parent',whImg_ax, 'XData',[0 size(whImg,2)*pxSzT]);
    if isempty(cropROIpos)
        cropROIpos = [whImg_ax.XLim(1) whImg_ax.YLim(1) ...
                      whImg_ax.XLim(2) whImg_ax.YLim(2)];
    end
    rectangle('Position', cropROIpos, 'Parent',whImg_ax,...
        'EdgeColor','w', 'LineWidth',3)
    colormap(parula(256))
    set(whImg_ax, 'XTick',[], 'FontSize',whImg_axFontSz)
    set(whImg_ax, 'XLim',[0 ...
        max(size(whImg,2)*pxSzT+delayImgToCurr,xt_cur(end))-delayImgToCurr])
    title('whole image, filtered')
    %xlabel('t (ms)')
    ylabel('x (pixels)')
 
else
    whImg_ax = axes('Parent',whImg_fig_final, ...
        'Units','normalized', ...
        'Position',[0.03 0.55 0.94 0.22]);
    image(whImg, 'CDataMapping','scaled', ...
        'Parent',whImg_ax, 'XData',[0 size(whImg,2)*pxSzT]);
    if isempty(cropROIpos)
        cropROIpos = [whImg_ax.XLim(1) whImg_ax.YLim(1) ...
            whImg_ax.XLim(2) whImg_ax.YLim(2)];
    end
    rectangle('Position',cropROIpos, 'Parent',whImg_ax,...
        'EdgeColor','w', 'LineWidth',5)
    colormap(whImg_ax, parula(256))
    set(whImg_ax, 'XTick',[], 'FontSize',whImg_axFontSz)
    title('whole image, filtered')
    %xlabel('t (ms)')
    ylabel('x (pixels)')
    delayImgToCurr = 0;
    xt_cur = 0;
    % colorbar   
    clim(whImg_ax, getAxisLimits(whImg, 1))
    h_cb_whImg = colorbar(whImg_ax,'north');
    h_cb_whImg.Position = [0.90 0.78 0.07 0.01];
    h_cb_whImg.Label.String = 'F';
end

whImgProf = axes('Parent',whImg_fig_final, ...
    'Units','normalized', ...
    'Position',[0.03 0.415 0.94 0.13]);
plot(linspace(0, size(whImg,2)*pxSzT, size(whImg,2)), mean(whImg,1), ...
    'Parent',whImgProf, 'Color','b');
xlabel(whImgProf,'t (ms)')
ylabel(whImgProf,'average fluorescence')
set(whImgProf,'FontSize',whImg_axFontSz)
set(whImgProf, 'XLim',[0 ...
    max(size(whImg,2)*pxSzT+delayImgToCurr,xt_cur(end))-delayImgToCurr],...
    'YLim',getAxisLimits(mean(whImg,1),5))

cImg_ax = axes('Parent',whImg_fig_final, ...
    'Units','normalized', ...
    'Position',[0.03 0.15 0.94 0.2]);
image(crImg, 'CDataMapping','scaled', 'Parent',cImg_ax, ...
    'YData',[0 size(crImg,1)*pxSzX], ...
    'XData',[0 size(crImg,2)*pxSzT]);
colormap(parula(256))
set(cImg_ax,'FontSize',whImg_axFontSz) 
xlabel('t (ms)')
ylabel('x (\mum)')
title('cropped image (white rectangle), filtered & normalized')

% colorbar
% cLims = [floor(prctile(crImg(:),1)*10)/10 ceil( prctile(crImg(:),99.9) )];
clim(cImg_ax, getAxisLimits(prctile(crImg(:), [1,99.5]), 1))
% clim(cImg_ax, getAxisLimits(crImg, 1))
h_cb = colorbar(cImg_ax,'south');
h_cb.Position = [0.90 0.11 0.07 0.01];
h_cb.Label.String = '\DeltaF/F_0';
%h_cb.Label.Position = [1 2.8 0];
% show ROIs for selected profiles
for i = 1:height(selectedROIs)      
    hPatch = selectedROIs.patch(i);
    rectYData = hPatch.YData .* pxSzX;
    rectXData = hPatch.XData;
    rectangle('Position',[rectXData(1) rectYData(3) rectXData(2)-rectXData(1) rectYData(1)-rectYData(3)],'Parent',cImg_ax,...
          'EdgeColor','k','LineWidth',1)
    %copyobj(hPatch,cImg_ax);    
end 

%invisible axes for text 
text_ax = axes('Parent',whImg_fig_final, ...
    'Units','normalized', ...
    'Position',[0.03 0.01 0.94 0.11]);
set(text_ax,'Visible','off')

notes = hObjs.h_table_notes.Data;
lineNotes = strjoin(cellfun(@(x) [x,' // '], notes,'UniformOutput',0));

ImgDataPath = fullfile(imgData.filePath,imgData.fileName);
ImgDataPath = strrep(ImgDataPath, '_', '\_');

if isfield(getappdata(mainFig),'electroPhys')
    
    electroPhys = getappdata(mainFig,'electroPhys');
    
    if isfield(electroPhys,'cellParameters')
        
        cellParams = electroPhys.cellParameters;
        cellParamsStr = strjoin(cellfun(@(a,b,c) sprintf('%s = %0.2f %s; ',a,b,c),...
            cellParams.Properties.VariableNames,...
            num2cell(cellParams{1,:}),...
            cellParams.Properties.VariableUnits,...
            'UniformOutput',0));
        
        txtStr = {sprintf('animal: %s',hObjs.popUpMenuAnimal.String{hObjs.popUpMenuAnimal.Value});...
                  sprintf('%s',lineNotes);...
                  ['cell parameters: ',cellParamsStr];...
                  sprintf('average ICa = %0.2f (pA); average current density = %0.2f (pA/pF)',electro{end,3},electro{end,4});...
                  sprintf('ImgDataPath: %s',ImgDataPath)};
    else
        txtStr = {sprintf('animal: %s',hObjs.popUpMenuAnimal.String{hObjs.popUpMenuAnimal.Value});...
                  sprintf('%s',lineNotes);...
                  sprintf('average ICa = %0.2f (pA)',electro{end,3});...
                  sprintf('ImgDataPath: %s',ImgDataPath)};
    end
    
else
    txtStr = {sprintf('animal: %s',hObjs.popUpMenuAnimal.String{hObjs.popUpMenuAnimal.Value});...
              sprintf('%s',lineNotes);...
              sprintf('ImgDataPath: %s',ImgDataPath)};
       
end

text(min(get(text_ax,'XLim')), max(get(text_ax,'YLim')), txtStr,...
    'Parent',text_ax, 'FontUnits','points', 'FontSize',14, ...
    'VerticalAlignment','top')

% figure of xy image, position of scanning line and TPP point, scale 
if ~isempty(imgData.imgDataXY)
    scanLinePos = imgData.scanLinePos;
    xyCellImg = imgData.imgDataXY;
    
    rp = min(scRes(3:4))/max(scRes(3:4));    
    dres = min(scRes(3:4))*2 - max(scRes(3:4));   
    szAx = min(scRes(3:4))-floor(dres/2);
    % 1.ch      
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
    line([ch1_h.XLim(2)-5-10/pxSzX ch1_h.XLim(2)-5], ...
        [ch1_h.YLim(2)-5 ch1_h.YLim(2)-5], ...
        'Parent',ch1_h, 'LineWidth',2, 'Color','g', 'LineStyle','-'); 
    text(ch1_h.XLim(2)-5-5/pxSzX, ch1_h.YLim(2)-15, ...
        '10 \mum', 'HorizontalAlignment','center',...
        'Color','g', 'FontSize',12, 'Parent',ch1_h, 'FontWeight','bold')
    
    % 2.ch
    try   
    ch2_h = axes('Parent',xyCellImg_h, 'Units','normalized', ...
        'Position',[0.6/(scRes(3)/scRes(4)) + 0.1 0.2 0.6/(scRes(3)/scRes(4)) 0.6]);
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

% get handles of figures of spark recovery analysis
h_figsSparkRecAnalysis = ...
    findobj('-regexp','Name','repetitive sparks analysis');

% handles of figures of ca events  
h_figsCaEvents = findobj('Type','figure','Tag','CaEventParam');
if ~isempty(h_figsCaEvents)
    eventNum = arrayfun(@(x) sscanf(x.Name,'posOfRoi %*f \\mum; spark #%d'),...
        h_figsCaEvents,'UniformOutput',1);
    posOfROI = arrayfun(@(x) sscanf(x.Name,'posOfRoi %f \\mum; spark #%*d'),...
        h_figsCaEvents,'UniformOutput',1);
    [~,indx_e] = sortrows([posOfROI,eventNum],[1 2]);
    h_figsCaEvents = h_figsCaEvents(indx_e);
end

% get handles of all figures to save
if hObjs.check_saveEventsFigs.Value
    allHtoSave = [whImg_fig_final; 
        h_figsSparkRecAnalysis; xyCellImg_h; h_figsCaEvents];
else
    allHtoSave = [whImg_fig_final; h_figsSparkRecAnalysis; xyCellImg_h];
    close(h_figsCaEvents)
end

% get path to save figures 
[pathFigs, nameFigs, ~] = fileparts(path_xls);

% delete files if exist
if exist(sprintf('%s/%s.pdf',pathFigs,nameFigs), 'file')
    delete(sprintf('%s/%s.pdf',pathFigs,nameFigs))
end

if exist(sprintf('%s/%s.ps',pathFigs,nameFigs), 'file')
    delete(sprintf('%s/%s.ps',pathFigs,nameFigs))
end

% save figures as .pdf
if ~isempty(allHtoSave)
    for i = 1:numel(allHtoSave)
              
        % % set upfigures for printing
        set(allHtoSave(i),'Units','pixels');
        posFig = get(allHtoSave(i),'Position');
        set(allHtoSave(i),'PaperPositionMode','Auto',...
                          'PaperUnits','points',...
                          'PaperSize',[posFig(3), posFig(4)])
        % 
       % create file 
        if ~exist(sprintf('%s/%s.pdf',pathFigs,nameFigs), 'file')
            % print(allHtoSave(i), ...
            %     sprintf('%s/%sAllEvents.ps',pathFigs,nameFigs), '-dpsc')
            exportgraphics(allHtoSave(i), ...
                sprintf('%s/%s.pdf',pathFigs,nameFigs))
        else
            % print(allHtoSave(i), '-dpsc', '-append', ...
            %     sprintf('%s/%sAllEvents.ps',pathFigs,nameFigs))
            exportgraphics(allHtoSave(i), ...
                sprintf('%s/%s.pdf',pathFigs,nameFigs), ...
                'Append',true)
        end                             
    end   
    close(allHtoSave)
end

%%
close(findobj('Type','figure','Name','Experiment prewiev'));

% clear window
cla(hObjs.ax_img_sparks);
cla(hObjs.ax_img);
cla(hObjs.h_ax_transCh);

delete(get(hObjs.ax_img,'Children'))
cla(hObjs.ax_prof);

set(mainFig,'Pointer','arrow')
drawnow

% if this re-analyzed image, called from spark recovery analysis than close
% window
if isfield(getappdata(mainFig),'reAnalyzeImg')
    if getappdata(mainFig,'reAnalyzeImg') == true
        close(mainFig)
    end
end

% % remove after
% reAnalyze([],[],mainFig)

end






