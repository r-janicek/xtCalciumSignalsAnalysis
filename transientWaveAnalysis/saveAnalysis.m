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

filter = fullfile(char(imgData.filePath),char([name,'_analysisResults','.xls']));      
[file,path] = uiputfile(filter);

if path == 0
    return
end 

set(mainFig,'Pointer','watch')
drawnow

% set up path for xls file
path_xls = fullfile(path,file);


% % %%%%%%%%%%
% keyboard
% % % save figure for pulication
%  dataToWrite = imgData.imgDataXTfluoFN;
% 
%  dataToWriteF = wiener2(dataToWrite,[3 3]);
%  dataToWriteF = imgaussfilt(dataToWriteF,0.5);
%  % no negative numbers
%  dataToWriteF = dataToWriteF + abs(min(dataToWriteF(:)));
% % 
% % 
% hf = figure('units','normalized','outerposition',[0 0.05 1 0.95]);
% axI = axes('Parent',hf);
% set(axI,'Position',hObjs.ax_img.Position)
% imagesc(axI,dataToWriteF)
% colormap(axI,parula(512))
% 
% axI.CLim = [0,8];
% caxis(axI,[0,8])
% hcb = colorbar(axI,'northoutside');
% 
% % copy profile
% copyobj(hObjs.ax_prof,hf)
% 
% % copy scale in um
% copyobj(hObjs.ax_sc,hf)
% 
% % set upfigures for printing
% set(hf,'Units','pixels');
% posFig = get(hf,'Position');
% set(hf,'PaperPositionMode','Auto',...
%     'PaperUnits','points',...
%     'PaperSize',[posFig(3), posFig(4)])
% 
% print(hf,fullfile(path,[name,'_cropped']),'-dpdf','-r600')
% 
% %%%%%%%%%%


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
    calcParamsFromFit(2:end,[1,3,4,5,13]) = ...
        num2cell(cell2mat(calcParamsFromFit(2:end,[1,3,4,5,13])) + posOfROIs(1));
    
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
posOfCropROI = [{'position of cropROI'};{'x(ms); y(px); w(ms); h(px)'};...
    num2cell(imgData.cropROIpos)'];

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

TPP_intervals = imgData.TPP_delays;  
               
TPP_delays = {'intra-pair interval of photolytic pulses (ms)', TPP_intervals(1,1);...
              'pair-to-pair interval of photolytic pulses (ms)',TPP_intervals(2,1)};
durTPP = {'duration of photolitic pulse (ms)',imgData.durOfTPP};



%% save electrophysiology data, only when data preview is used

if isfield(getappdata(mainFig),'electroPhys')
    
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
     
end


%% save data to xls file
% add folder where java classes for creating of xls files are
currentFolder = pwd;
expression = 'Matlab';
splitStr = regexp(currentFolder,expression,'split');

javaaddpath([splitStr{1,1},'Matlab/createXLSonMAC/poi_library/poi-3.8-20120326.jar']);
javaaddpath([splitStr{1,1},'Matlab/createXLSonMAC/poi_library/poi-ooxml-3.8-20120326.jar']);
javaaddpath([splitStr{1,1},'Matlab/createXLSonMAC/poi_library/poi-ooxml-schemas-3.8-20120326.jar']);
javaaddpath([splitStr{1,1},'Matlab/createXLSonMAC/poi_library/xmlbeans-2.3.0.jar']);
javaaddpath([splitStr{1,1},'Matlab/createXLSonMAC/poi_library/dom4j-1.6.1.jar']);
javaaddpath([splitStr{1,1},'Matlab/createXLSonMAC/poi_library/stax-api-1.0.1.jar']);

% delete xls file if exist
if exist(path_xls, 'file')
    delete(path_xls)
end

% save img and analysis info 
xlwrite(path_xls, imgPath, 'info', 'A1');
xlwrite(path_xls, imgName, 'info', 'A2');

xlwrite(path_xls, imgSize_px, 'info', 'A4');
xlwrite(path_xls, imgSize, 'info', 'A6');
xlwrite(path_xls, pxSz_x, 'info', 'A8');
xlwrite(path_xls, pxSz_t, 'info', 'A9');
xlwrite(path_xls, blank, 'info', 'A10');

xlwrite(path_xls, durTPP, 'info', 'A13');
xlwrite(path_xls, TPP_delays, 'info', 'A14');

%%%%%%%
if isfield(getappdata(mainFig),'electroPhys')
    xlwrite(path_xls, delayImg, 'info', 'A17');
end
%%%%%%
xlwrite(path_xls, animal, 'info', 'A19');
xlwrite(path_xls, exp_notes, 'info', 'A20');

xlwrite(path_xls, imgFiltersUsed, 'info', 'C18');

xlwrite(path_xls, posOfCropROI, 'info', 'C23');

xlwrite(path_xls, result, 'info', 'E1');


% check whether save profiles and fits
if hObjs.check_saveProfsAndFits.Value
    
    % save all profiles data
    for i=1:numel(profiles)
        nCol = size(profiles{i},2);
        xlwrite(path_xls, profiles{i}, 'profilesOfEvents', [xlscol((i-1)*nCol+1),'1']);
    end
                 
end

% save current ...
if exist('electro','var')
    xlwrite(path_xls, electro, 'electroPhys', 'A1');
end

if exist('cellParameters','var')
    
    xlwrite(path_xls,[cellParameters.Properties.VariableNames;...
        cellParameters.Properties.VariableUnits;...
        table2cell(cellParameters)], 'currentAndCell', 'G1');
    
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
d_TP_laser = imgData.d_TP_laser;

% plot whole image with white rectangle where it was cropped
if isfield(getappdata(mainFig),'electroPhys') 
    
    electroPhys = getappdata(mainFig,'electroPhys');
    
    IcaLimits = electroPhys.IcaLimits;   
    delayImgToCurr = electroPhys.delayImgToCurr; % in ms
    xt_cur = electroPhys.xt_cur;
    V = electroPhys.voltage;
    ICa = electroPhys.current;
    [~,pos_delay] = min(abs(xt_cur-delayImgToCurr));
       
    vol_ax = axes('Parent',whImg_fig_final,'Units','normalized','Position',[0.04 0.93 0.94 0.06]);
    plot(xt_cur(pos_delay:end),V(pos_delay:end),'Parent',vol_ax,'Color','k','LineWidth',2);
    ylabel(vol_ax,'voltage (mV)')
    set(vol_ax,'XTick',[],'FontSize',14)
    set(vol_ax,'XLim',[xt_cur(pos_delay) max(size(whImg,2)*pxSzT+delayImgToCurr,xt_cur(end))],...
        'YLim',[min(V(pos_delay:end))-5 max(V(pos_delay:end))+5])
    
    curr_ax = axes('Parent',whImg_fig_final,'Units','normalized','Position',[0.04 0.795 0.94 0.13]);
    plot(xt_cur(pos_delay:end),ICa(pos_delay:end),'Parent',curr_ax,'Color','r');
    ylabel(curr_ax,'current (pA)')
    set(curr_ax,'XTick',[],'FontSize',14)
    set(curr_ax,'XLim',[xt_cur(pos_delay) max(size(whImg,2)*pxSzT+delayImgToCurr,xt_cur(end))],...
        'YLim',IcaLimits)
    
    whImg_ax = axes('Parent',whImg_fig_final,'Units','normalized','Position',[0.04 0.55 0.94 0.22]);
    image(whImg,'CDataMapping','scaled','Parent',whImg_ax,'XData',[0 size(whImg,2)*pxSzT]);
    if isempty(cropROIpos)
        cropROIpos = [whImg_ax.XLim(1) whImg_ax.YLim(1) whImg_ax.XLim(2) whImg_ax.YLim(2)];
    end
    rectangle('Position',cropROIpos,'Parent',whImg_ax,...
        'EdgeColor','w','LineWidth',3)
    colormap(parula(256))
    set(whImg_ax,'XTick',[],'FontSize',14)
    set(whImg_ax,'XLim',[0 max(size(whImg,2)*pxSzT+delayImgToCurr,xt_cur(end))-delayImgToCurr])
    title('whole image, filtered')
    %xlabel('t (ms)')
    ylabel('x (pixels)')
 
else
    whImg_ax = axes('Parent',whImg_fig_final,'Units','normalized','Position',[0.04 0.55 0.94 0.22]);
    image(whImg,'CDataMapping','scaled','Parent',whImg_ax,'XData',[0 size(whImg,2)*pxSzT]);
    if isempty(cropROIpos)
        cropROIpos = [whImg_ax.XLim(1) whImg_ax.YLim(1) whImg_ax.XLim(2) whImg_ax.YLim(2)];
    end
    rectangle('Position',cropROIpos,'Parent',whImg_ax,...
        'EdgeColor','w','LineWidth',5)
    colormap(parula(256))
    set(whImg_ax,'XTick',[],'FontSize',14)
    title('whole image, filtered')
    %xlabel('t (ms)')
    ylabel('x (pixels)')
    delayImgToCurr = 0;
    xt_cur = 0;
    
end

% copy filtered and normalized cropped image
ax_CroppedImg = copyobj(hObjs.ax_img,whImg_fig_final);
ax_CroppedImg.Position = [0.04 0.305 0.94 0.22];
title(ax_CroppedImg,'cropped image (white rectangle), filtered & normalized')
% add colorbar
colormap(ax_CroppedImg,parula(256))
h_CroppedImg = findall(ax_CroppedImg,'Type','Image');
caxis(ax_CroppedImg,[ floor(prctile(h_CroppedImg.CData(:),1)*10)/10 ...
    ceil(prctile(h_CroppedImg.CData(:),99.9)) ])
h_cb = colorbar(ax_CroppedImg,'south');
h_cb.Position = [0.92 0.035 0.06 0.01];
h_cb.Label.String = '\DeltaF/F_0';


% copy profile
ax_imgProf = copyobj(hObjs.ax_prof,whImg_fig_final);
ax_imgProf.Position = [0.04 0.11 0.94 0.19];

% show time to first wave from beginning of last triggered transient
m_triggTrans = strcmp(result(:,6), 'transient_trigg' );
m_wave = strcmp(result(:,6), 'wave' );

transRes = cell2mat(result(2:end,9));
transRes = [nan;transRes];
transRes(~m_triggTrans) = nan;

[lastTrans_t0,indTr]= max(transRes);
firstWave_t0 = min([result{m_wave,9}]);

try
    yL = result{indTr,10} + result{indTr,15};
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

%plot position of TPP
try   
    if isempty(crop_s_t)
        crop_s_t = 0;
    end
    
    s_TPP_whImg = imgData.s_TPP + crop_s_t;
    e_TPP_whImg = imgData.e_TPP + crop_s_t;
    s_TPP_c = s_TPP_whImg-crop_s_t;
    e_TPP_c = e_TPP_whImg-crop_s_t;
    durOfTPP = imgData.durOfTPP;
    posOfTPPinScanLine = imgData.posOfTPPinScanLine; % in whole image

    if ~isempty(s_TPP_whImg)
        
        if ~isempty(posOfTPPinScanLine)
            
            for i=1:numel(s_TPP_whImg)
                
                pos_c = [t(s_TPP_c(i)), (posOfTPPinScanLine - (d_TP_laser/pxSzX)/2 - cropROIpos(2))*pxSzX, durOfTPP, d_TP_laser];
                rectangle('Position',pos_c,'Parent',cImg_ax,'EdgeColor','r','LineWidth',2)
                
                pos = [t_whImg(s_TPP_whImg(i)), posOfTPPinScanLine - (d_TP_laser/pxSzX)/2, durOfTPP, d_TP_laser/pxSzX];
                rectangle('Position',pos,'Parent',whImg_ax,'EdgeColor','r','LineWidth',2)
                
            end
            
        else
            for i=1:numel(s_TPP_whImg)
                
                
                line([t(s_TPP_c(i)) t(s_TPP_c(i))],get(cImg_ax,'YLim'),'Parent',cImg_ax,'LineWidth',2,'Color','r',...
                    'LineStyle',':');
                line([t(e_TPP_c(i)) t(e_TPP_c(i))],get(cImg_ax,'YLim'),'Parent',cImg_ax,'LineWidth',2,'Color','r',...
                    'LineStyle',':');
                
                line([t_whImg(s_TPP_whImg(i)) t_whImg(s_TPP_whImg(i))],get(whImg_ax,'YLim'),'Parent',whImg_ax,'LineWidth',2,'Color','r',...
                    'LineStyle',':');
                line([t_whImg(e_TPP_whImg(i)) t_whImg(e_TPP_whImg(i))],get(whImg_ax,'YLim'),'Parent',whImg_ax,'LineWidth',2,'Color','r',...
                    'LineStyle',':');
            end
        end
    end
catch
end


% invisible axes for text 
text_ax = axes('Parent',whImg_fig_final,'Units','normalized','Position',[0.03 0.01 0.94 0.08]);
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

text(min(get(text_ax,'XLim')),max(get(text_ax,'YLim')),txtStr,...
    'Parent',text_ax,'FontUnits','points','FontSize',16,...
    'VerticalAlignment','top')



%% figure of xy image, position of scanning line and TPP point, scale 
if ~isempty(imgData.imgDataXY)
    
    TPPpointPos = imgData.TPPpointPos; 
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
   
    if ~isempty(TPPpointPos)       
        pos = [TPPpointPos(1)-(d_TP_laser/pxSzX)/2 TPPpointPos(2)-(d_TP_laser/pxSzX)/2 (d_TP_laser/pxSzX) (d_TP_laser/pxSzX)];
        rectangle('Position',pos,'Curvature',[1 1],'Parent',ch1_h,'EdgeColor','r','LineWidth',2)
    end
    
    % scale 10 um
    line([ch1_h.XLim(2)-5-10/pxSzX ch1_h.XLim(2)-5],[ch1_h.YLim(2)-5 ch1_h.YLim(2)-5],'Parent',ch1_h,'LineWidth',2,'Color','g',...
            'LineStyle','-'); 
    text(ch1_h.XLim(2)-5-5/pxSzX,ch1_h.YLim(2)-15,'10 \mum','HorizontalAlignment','center',...
        'Color','g','FontSize',12,'Parent',ch1_h,'FontWeight','bold')
    try   
    ch2_h =  axes('Parent',xyCellImg_h,'Units','normalized','Position',[0.6/(scRes(3)/scRes(4)) + 0.1 0.2 0.6/(scRes(3)/scRes(4)) 0.6]);
    image((xyCellImg{2,1}),'Parent',ch2_h,'CDataMapping','scaled')
    colormap(ch2_h,'gray')
    set(ch2_h,'XTick',[],'YTick',[])
    if ~isempty(scanLinePos)
        line([scanLinePos(1) scanLinePos(2)],[scanLinePos(3) scanLinePos(4)],'Parent',ch2_h,'LineWidth',2,'Color','y',...
            'LineStyle','-');        
    end
    if ~isempty(TPPpointPos)
        pos = [TPPpointPos(1)-(d_TP_laser/pxSzX)/2 TPPpointPos(2)-(d_TP_laser/pxSzX)/2 (d_TP_laser/pxSzX) (d_TP_laser/pxSzX)];
        rectangle('Position',pos,'Curvature',[1 1],'Parent',ch2_h,'EdgeColor','r','LineWidth',2)     
    end
    % scale 5 um
    line([ch1_h.XLim(2)-5-10/pxSzX ch1_h.XLim(2)-5],[ch1_h.YLim(2)-5 ch1_h.YLim(2)-5],'Parent',ch2_h,'LineWidth',2,'Color','g',...
            'LineStyle','-'); 
    text(ch1_h.XLim(2)-5-5/pxSzX,ch1_h.YLim(2)-15,'10 \mum','HorizontalAlignment','center',...
        'Color','g','FontSize',12,'Parent',ch2_h,'FontWeight','bold')
    catch
    end
    
else
    xyCellImg_h = [];
end


%% save figures

% get handles of figures of individial events analysis
h_figsEventsAnalysis = findall(0,'Type','figure','Tag','imgsOfEventsFromAnalysis');
% make them visible again
arrayfun(@(x) set(x,'Visible','on'), h_figsEventsAnalysis)

% get handles of all figures to save
if hObjs.check_saveEventsFigs.Value
    if isempty(h_figsEventsAnalysis)
        allHtoSave = [whImg_fig_final;xyCellImg_h];
    else
        allHtoSave = [whImg_fig_final;h_figsEventsAnalysis;xyCellImg_h];
    end
        
else   
    allHtoSave = [whImg_fig_final;xyCellImg_h];
end
           
% get path to save figures 
[pathFigs,nameFigs,~] = fileparts(path_xls);

% delete files if exist
if exist(sprintf('%s/%sAllEvents.pdf',pathFigs,nameFigs), 'file')
    delete(sprintf('%s/%sAllEvents.pdf',pathFigs,nameFigs))
end

if exist(sprintf('%s/%sAllEvents.ps',pathFigs,nameFigs), 'file')
    delete(sprintf('%s/%sAllEvents.ps',pathFigs,nameFigs))
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
        if ~exist(sprintf('%s/%sAllEvents.ps',pathFigs,nameFigs), 'file')  
            print(allHtoSave(i), ...
                sprintf('%s/%sAllEvents',pathFigs,nameFigs), '-dpsc')
            exportgraphics(allHtoSave(i), ...
                sprintf('%s/%sAllEvents.pdf',pathFigs,nameFigs))
        else
            print(allHtoSave(i), '-dpsc', '-append', ...
                sprintf('%s/%sAllEvents',pathFigs,nameFigs))
            exportgraphics(allHtoSave(i), ...
                sprintf('%s/%sAllEvents.pdf',pathFigs,nameFigs), ...
                'Append',true, 'Resolution',300)
        end
        
                                    
        % export_fig(allHtoSave(i),sprintf('%s/%sAllEvents',pathFigs,nameFigs),'-pdf','-append')
                                     
    end
                
    close(allHtoSave)
end


% % create pdf (with resolution 300 dpi) from postscript file (need ghostscript)
% cmdToCreatePDF = sprintf( '/usr/local/bin/ps2pdf -r300 %s %s',...
%     fullfile(pathFigs,sprintf('%sAllEvents.ps',nameFigs)),...
%     fullfile(pathFigs,sprintf('%sAllEvents.pdf',nameFigs)) );
% 
% [status,cmdout] = system(cmdToCreatePDF);

% create pdf from postscript file (need ghostscript)
% ps2pdf('psfile', fullfile(sprintf('%s/%sAllEvents.ps',pathFigs,nameFigs)),...
%     'pdffile', fullfile(sprintf('%s/%sAllEvents.pdf',pathFigs,nameFigs)),...
%     'gspapersize', 'a4', 'verbose', 0,'resolution',300, ...
%             'gscommand', [splitStr{1,1},'Matlab/ghostscript/9.26/bin/gs'], ...
%             'gsfontpath', [splitStr{1,1},'Matlab/ghostscript/9.26/lib'], ...
%             'gslibpath', [splitStr{1,1},'Matlab/ghostscript/9.26/lib'])


%%
close(findobj('Type','figure','Name','analysis of selected ROIs'));

% clear window
cla(hObjs.ax_img);
cla(hObjs.ax_prof);

delete(get(hObjs.ax_img,'Children'))
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

