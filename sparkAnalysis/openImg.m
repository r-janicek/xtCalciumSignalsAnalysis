function openImg(hO, ~,mainFig)
% open image file using bioformats

% check how openImg function is called
if ~isempty(hO)
    if contains(hO.String,'open')
        % set batch processing mode
        setappdata(mainFig,'batchProcessing',false)
    end
end

% reAnalyze
if isfield(getappdata(mainFig),'reAnalyzeImg')
    if getappdata(mainFig,'reAnalyzeImg')==1
        reAnalyzeImgFlag = true;
    else
        reAnalyzeImgFlag = false;
    end
else    
    reAnalyzeImgFlag = false;    
end

% batch processing
if isfield(getappdata(mainFig),'batchProcessing')
    if getappdata(mainFig,'batchProcessing')
        batchProcessing = true;
    else
        batchProcessing = false;
    end
else    
    batchProcessing = false;    
end

% get path to file
if (isfield(getappdata(mainFig),'imgData') && reAnalyzeImgFlag) || batchProcessing % getappdata(mainFig,'loadFlag')==1

    filePath = getappdata(mainFig,'FilePath');
    fileName = getappdata(mainFig,'FileName');
    
    if ~iscell(fileName)
        fileName = {fileName};
          
    end
    
else
    
    if isfield(getappdata(mainFig),'imgData')
        
        imgDataOld = getappdata(mainFig,'imgData');
        pathOld = imgDataOld.filePath;
              
        [fileName, filePath] = uigetfile({'*.oif;*.oib;*.tif;*.pic'},'Load image',...
            pathOld,'Multiselect','off');
    else
        
        [fileName, filePath] = uigetfile({'*.oif;*.oib;*.tif;*.pic'},'Load image',...
            '/Users/radoslavjanicek/Desktop/',...
            'Multiselect','off');
    end
    
    if ~iscell(fileName)
        fileName = {fileName};
        
    end
    
end

if filePath == 0
   clear FilePath 
   return 
end

set(mainFig,'Pointer','watch')
drawnow

if isfield(getappdata(mainFig),'imgData') && isfield(getappdata(mainFig),'seriesOfImgs')
    rmappdata(mainFig,'imgData'); 
    rmappdata(mainFig,'seriesOfImgs'); 
end

[~,~,ext] = fileparts(fileName{1,1});

% load data using OME_bioformats
try
    % load all data about image using OME_bioformats
    data = bfopen(fullfile(filePath,fileName{1,1}));
 
    % meta data
    metaData = cat(2,cell(data{1,2}.keySet.toArray),cell(data{1,2}.values.toArray));
    
    indx_SM = find(~cellfun(@isempty,regexp(metaData(:,1),'\<ScanMod+e','match'))==1,1,'first');
    
    try
        scanMode = metaData{indx_SM,2};
    catch
        scanMode = [];
    end
    try
        for i = 1:size(metaData,1)
            
            ImageDescription_Names{i,1} = genvarname(char(metaData{i,1}));
            ImageDescription_val {i,1} = metaData{i,2};
            
        end
        
        [~,ia,~] = unique(ImageDescription_Names);
        metaDataS = cell2struct(cell(ImageDescription_val(ia)), cell(ImageDescription_Names(ia)),1);
    catch
        
    end
    %OME meta data
    omeMetaData = data{1,4};
    omeMetaDataXML = char(omeMetaData.dumpXML());
    
    omeMD = textscan(omeMetaDataXML,'%s','Delimiter',{'<','>'});
    omeMD1 = omeMD{1,1}(~cellfun(@isempty,regexp(omeMD{1,1},'=')));

    try
        microscopeVer = metaDataS.Global0x5BVersionInfo0x5DSystemName;
        userName = metaDataS.Global0x5BFileInfo0x5DUserName;
    catch 
        microscopeVer = '';
        userName = '';
    end
    
    % image data, change to double and rotate
    % chNames = arrayfun(@(x) sprintf('ch%d',x) ,(1:numel(data{1,1}(:,1))), 'UniformOutput',0);
    % imgDataXT = cell2struct(data{1,1}(:,1), chNames, 1);
    imgDataXT = data{1,1}(:,1);
    imgDataXT = cellfun(@(x) double(imrotate(x,90)),imgDataXT, 'UniformOutput',0);
    
    % assign channels and set up popup menus
    hObjs = getappdata(mainFig,'hObjs');
    chNames = arrayfun(@(x) sprintf('ch #%d',x) ,(1:numel(data{1,1}(:,1))), 'UniformOutput',0);
    chNames = [chNames,{'none'}];
    
    set(hObjs.popUpMenuChFluo,'String',chNames,'Value',1)
    set(hObjs.popUpMenuChTrans,'String',chNames,'Value',2)
         
    imgDataXTfluo = imgDataXT{1,1};
    
    try
        imgDataXTtrans = imgDataXT{2,1};
    catch
        imgDataXTtrans = [];
    end
    
    % save XY image data, if possible
    try
        imgDataXY = data{2,1}(:,1);
    catch
        imgDataXY = [];
    end
    
    % load pixel size, time dimension
    switch ext
        
        case '.pic'
            pxSzT = 1000/str2double(metaDataS.GlobalINFO_RASTER_RATE);
            
        case '.oib'
            pxSzT = str2double(metaDataS.GlobalTimePerLine)/1000;
            
        case '.oif'
            pxSzT = str2double(metaDataS.GlobalTimePerLine)/1000;
            
        case '.tiff'
            pxSzT = double(omeMetaData.getPlaneDeltaT(0,0).value);
            
    end
    
    pxSzX = double(omeMetaData.getPixelsPhysicalSizeX(0).value);
    
    if pxSzX>100000 % looks like error in olympus writing metadata
        pxSzX = pxSzX/1000000;
    end
    
    t = linspace(0,(size(imgDataXTfluo,2)-1)*pxSzT,size(imgDataXTfluo,2)); %ms
    
    % try to load positions of ROIs from image metadata
    if strcmp(ext,'.oif')
        
        fileList = dir([fullfile(filePath,fileName{1,1}),'.files']);
        fNnames = {fileList.name};
        roiFile = fNnames{cellfun(@(x) ~isempty(x),strfind(fNnames,'.roi'))==1};
        roiFile = fullfile([fullfile(filePath,fileName{1,1}),'.files'],roiFile);
        
        fileID = fopen(roiFile);
        fileROIlines = textscan(fileID,'%s','Delimiter','\n');
        fileROIlines = regexprep(fileROIlines{1,1},'[^a-zA-Z_-.,=0-9]','');
        fclose(fileID);
        
        numOfROIs = find(cellfun(@(x) ~isempty(x),strfind(fileROIlines,'ROIBaseFileInformation'))==1);
        
        for i=1:length(numOfROIs)
            
            try
                ROI = fileROIlines(numOfROIs(i):numOfROIs(i+1)-1);
            catch
                ROI = fileROIlines(numOfROIs(i):end);
            end
            
            posOfLine = find((~cellfun(@isempty,regexpi(ROI,'SHAPE'))));
            shapeNum = regexpi(ROI{posOfLine,1},'\d*','match');
            
            switch char(shapeNum)
                
                case '3' % line
                    %keyboard
                    posOfX = (~cellfun(@isempty,regexp(ROI,'^X\w*')));
                    posOfY = (~cellfun(@isempty,regexp(ROI,'^Y\w*')));
                    
                    posOfX = sscanf(ROI{posOfX}, 'X=%f,%f');
                    posOfY = sscanf(ROI{posOfY}, 'Y=%f,%f');
                    
                    pos_l = [posOfX;posOfY];
                    pos_l = pos_l+1;
                    
                case '2' % point
                    
                    posOfX = (~cellfun(@isempty,regexp(ROI,'^X\w*')));
                    posOfY = (~cellfun(@isempty,regexp(ROI,'^Y\w*')));
                    
                    posOfX = sscanf(ROI{posOfX}, 'X=%f,%f');
                    posOfY = sscanf(ROI{posOfY}, 'Y=%f,%f');
                    
                    pos_P = [posOfX;posOfY];
                    pos_P = pos_P+1;
                    
            end
            
            clearvars ROI
            
        end
        
    else
        scanline = char(omeMD1(~cellfun(@isempty,regexp(omeMD1,'Line'))));
        pointTP = char(omeMD1(~cellfun(@isempty,regexp(omeMD1,'Point'))));
        
        scanline = char(regexp(scanline,'X1=".*" X2=".*" Y1=".*" Y2=".*"','match'));
        pointTP = char(regexp(scanline,'X=".*" Y=".*"','match'));
        
        pos_l = sscanf(scanline, 'X1="%f" X2="%f" Y1="%f" Y2="%f"');
        pos_P = sscanf(pointTP, 'X="%f" Y="%f"');
        pos_l = pos_l+1;
        pos_P = pos_P+1;
        
    end
    
    if ~exist('pos_l','var')
        pos_l = [];
    end
    
    if ~exist('pos_P','var')
        pos_P = [];
    end
    
    %filter data and plot image
    [imgDataXTfluo,ImgFiltersUsed] = imgFiltering(imgDataXTfluo,pxSzT,pxSzX);
    
    ax_img = hObjs.ax_img;
    cla(ax_img)
    
    image(imgDataXTfluo,'YData',[1 size(imgDataXTfluo,1)],'XData',[0 max(t)],...
        'CDataMapping','scaled','Parent',ax_img);
    set(ax_img,'XTick',[],'FontSize',14)
    set(get(ax_img,'Ylabel'),'String','x (px)','FontWeight','bold')
    
    % time profile
    prof_t = mean(imgDataXTfluo,1);
    
    % clear axes
    ax_img_sparks = hObjs.ax_img_sparks;
    cla(ax_img_sparks)
    % ax_img_sparks_2 = hObjs.ax_img_sparks_2;
    % cla(ax_img_sparks_2)
    
    ax_prof = hObjs.ax_prof;
    plot(t,prof_t,'Parent',ax_prof);
    set(ax_prof,'FontSize',14)
    set(get(ax_prof,'Xlabel'),'String','t (ms)','FontWeight','bold')
    set(get(ax_prof,'Ylabel'),'String','Fluorescence (F)','FontWeight','bold')
    xlim(ax_prof,[min(t) max(t)])
    
    try
        ylim(ax_prof,[min(prof_t) max(prof_t)])
    catch
        ylim(ax_prof,[min(prof_t) max(prof_t)+1])
    end
    
    % link axes
    linkaxes([ax_img,ax_img_sparks,ax_prof],'x')
    
    % load positions of photolytic pulses, from trans channel
    if strcmp(getappdata(mainFig,'analysisType'),'spark recovery photolysis')
        
        [s_TPP,e_TPP,durOfTPP,TPP_delays,...
            d_TP_laser,posOfTPPinScanLine] = loadPhotolysisPositions(imgDataXTfluo,...
                imgDataXTtrans,...
                t,...
                pxSzT,...
                pxSzX,...
                pos_P,...
                pos_l,...
                hObjs.ax_img,...
                1);
             
         % set as max of image in found positions of photolytic pulses
         % just to prevent, that when cropped it might miss position of
         % photolytic pulse, mainly when pulse is very short ~1ms and
         % scanning speed slow: time per line more than 1 ms
         for idx=1:numel(s_TPP)
            imgDataXTtrans(:,s_TPP(idx)+1:e_TPP(idx)) = max(imgDataXTtrans(:));
         end 
            
    else
        
        TPP_delays = [nan;nan];
        posOfTPPinScanLine = [];
        durOfTPP = [];
        s_TPP = [];
        e_TPP = [];
        d_TP_laser = 0;
    end
        
catch
    try
        data = imread(fullfile(filePath,fileName{1,1}));
        if ~iscell(data)
            data = {data};
        end
        metaData = imfinfo(fullfile(filePath,fileName{1,1}));
        fieldsMetaData = fieldnames(metaData);
        
        imgDataXT = data;
        imgDataXT = cellfun(@(x) double(imrotate(x,90)),imgDataXT, 'UniformOutput',0);

        imgDataXTfluo = imgDataXT{1,1};
        imgDataXY = [];
        imgDataXTtrans = [];
        
        TPP_delays = [nan;nan];
        posOfTPPinScanLine = [];
        durOfTPP = [];
        s_TPP = [];
        e_TPP = [];
        d_TP_laser = 0;
        pos_P = [];
        pos_l = [];
        pxSzT = 1;
        pxSzX = 1;
         
        t = linspace(0,(size(imgDataXTfluo,2)-1)*pxSzT,size(imgDataXTfluo,2)); %ms
        
        %filter data and plot image
        [imgDataXTfluo,ImgFiltersUsed] = imgFiltering(imgDataXTfluo,pxSzT,pxSzX);
        
        % assign channels and set up popup menus
        hObjs = getappdata(mainFig,'hObjs');
        chNames = arrayfun(@(x) sprintf('ch #%d',x) ,(1:numel(data)), 'UniformOutput',0);
        chNames = [chNames,{'none'}];
    
        set(hObjs.popUpMenuChFluo,'String',chNames,'Value',1)
        set(hObjs.popUpMenuChTrans,'String',chNames,'Value',2)
        
        % show image
        ax_img = hObjs.ax_img;
        cla(ax_img)
        
        image(imgDataXTfluo,'YData',[1 size(imgDataXTfluo,1)],'XData',[0 max(t)],...
            'CDataMapping','scaled','Parent',ax_img);
        set(ax_img,'XTick',[],'FontSize',14)
        set(get(ax_img,'Ylabel'),'String','x (px)','FontWeight','bold')
        
        % time profile
        prof_t = mean(imgDataXTfluo,1);
        
        % clear axes
        ax_img_sparks = hObjs.ax_img_sparks;
        cla(ax_img_sparks)
        % ax_img_sparks_2 = hObjs.ax_img_sparks_2;
        % cla(ax_img_sparks_2)
        
        ax_prof = hObjs.ax_prof;
        plot(t,prof_t,'Parent',ax_prof);
        set(ax_prof,'FontSize',14)
        set(get(ax_prof,'Xlabel'),'String','t (ms)','FontWeight','bold')
        set(get(ax_prof,'Ylabel'),'String','Fluorescence (F)','FontWeight','bold')
        xlim(ax_prof,[min(t) max(t)])
        
        try
            ylim(ax_prof,[min(prof_t) max(prof_t)])
        catch
            ylim(ax_prof,[min(prof_t) max(prof_t)+1])
        end
        
        % link axes
        linkaxes([ax_img,ax_img_sparks,ax_prof],'x')
        
        
    catch exception
        errordlg('can not open image file!')
    end
end


% maybe use percentile 1 and 99
z_max_img = max(imgDataXTfluo(:));
z_min_img = min(imgDataXTfluo(:));

% remove old image data
if isfield(getappdata(mainFig),'imgData')    
    rmappdata(mainFig,'imgData');    
end

% save image data
imgData = struct('filePath',filePath,...
                 'fileName',fileName,...
                 'imgDataXT',{imgDataXT},...
                 'imgDataXY',{imgDataXY},...
                 'wholeImgXT',{imgDataXT},...
                 'wholeImgFluoXT',imgDataXTfluo,...   % image data, filtered
                 'imgDataXTfluoR',imgDataXT{1,1},...  % raw image data
                 'imgDataXTfluoRN',imgDataXT{1,1},... % raw image data, normalized
                 'imgDataXTfluoF',imgDataXTfluo,...   % image data, filtered
                 'imgDataXTfluoFN',imgDataXTfluo,...   % image data, filtered and normalized
                 'imgDataXTtrans',imgDataXTtrans,... % image data transmitted light
                 'pxSzT',pxSzT,...
                 'pxSzX',pxSzX,...
                 't',t,...
                 'imgFiltersUsed',ImgFiltersUsed,...
                 'z_max_img',z_max_img,...
                 'z_min_img',z_min_img,...
                 's_TPP',s_TPP,...
                 'e_TPP',e_TPP,...
                 'durOfTPP',durOfTPP,...
                 'TPP_delays',TPP_delays,...
                 'd_TP_laser',d_TP_laser,...
                 'TPPpointPos',pos_P,...
                 'scanLinePos',pos_l,...
                 'posOfTPPinScanLine',posOfTPPinScanLine,...
                 'microscopeVer',microscopeVer,...
                 'userName',userName);

setappdata(mainFig,'imgData',imgData)

% set up main window
set(hObjs.txt_name_img,'String',sprintf('%s',fullfile(filePath,char(fileName))),...
    'Tooltip',sprintf('%s',fullfile(filePath,char(fileName))));

% show trans image in small window
image(imgDataXTtrans,'CDataMapping','scaled','Parent',hObjs.h_ax_transCh,...
    'ButtonDownFcn',{@axesTransChFunction});
set(hObjs.h_ax_transCh,'XTick',[],'YTick',[],'PickableParts','all')
set(hObjs.h_ax_transCh,'ButtonDownFcn',{@axesTransChFunction})

switch getappdata(mainFig,'analysisType')
    case 'spark detection'
        set(hObjs.txt_SpF,'String',['# sp*100',char(956),'m-1*s-1'])
        set(hObjs.txt_correctedSpF,'String',{'#';['sp*100',char(956),'m-1*s-1']})
        
    otherwise
        set(hObjs.h_table_profs,'Data',zeros(5,2))
        % set normalization panel
        set(hObjs.h_pb_norm,'Enable','on')
        set(hObjs.h_pb_norm,'String','normalize data')
        set(hObjs.h_pb_norm,'Callback',{@norm_data,mainFig})
        set(hObjs.h_pb_norm,'FontWeight','normal')
        
        % set up panel for spark recovery
        % set width of ROI for spark recovery to 5 um width in pixels (odd number)
        n_px = round(6/pxSzX);
        if mod(n_px,2)==0
            n_px = n_px+1;
        end
        
        set(hObjs.h_edit_ROI,'String',num2str(n_px))
        
        % set up panel for spark detection
        % set(hObjs.check_visRect,'Value',1)
        set(hObjs.txt_SpF,'String',['# sp*100',char(956),'m-1*s-1'])
        
end

% check for position of UVflash, flash lamp
if hObjs.h_check_UVflash.Value
    getUVflashPos([],[],mainFig);
end

% set blank panel
set(hObjs.h_push_BlankROI,'String','set ROI for blank','FontWeight','normal')
set(hObjs.h_push_BlankROI,'Callback', {@setBlankROI,mainFig})
set(hObjs.h_push_BlankROI,'Enable','on')

set(hObjs.h_edit_Blank,'String','190') % 'nan' 190 152)
set(hObjs.h_edit_Blank,'Enable','on')

if isempty(imgDataXY)    
    set(hObjs.h_pb_GetBlank,'Enable','off')    
else
    set(hObjs.h_pb_GetBlank,'Enable','on')   
end

% set up pushbutton for cropping
set(hObjs.h_pb_crop,'String','<html> <p align="center"> set ROI <br> for cropping <html>')
set(hObjs.h_pb_crop,'FontWeight','normal')

if strcmp(get(hObjs.h_pb_crop,'String'),'crop')  
    
    set(hObjs.h_pb_crop,'String','set ROI for cropping')
    set(hObjs.h_pb_crop,'FontWeight','normal')
    
end
% set FWHM for spark detection, in pixels, 2 um
set(hObjs.h_edit_fFWHM,'String','10')
% set(hObjs.h_edit_fFWHM,'String',num2str(ceil(2/pxSzX)))

% remove old analysis
if isfield(getappdata(mainFig),'lastPressedPusbutton')
    rmappdata(mainFig,'lastPressedPusbutton');   
end

if isfield(getappdata(mainFig),'profileAnalysis')
    rmappdata(mainFig,'profileAnalysis');   
end

if isfield(getappdata(mainFig),'sparkDetection')
    rmappdata(mainFig,'sparkDetection');   
end

% restore image slider
imgSld_w = str2double(hObjs.h_edit_imgSld_w.String)*1000; % change from s to ms
set(hObjs.h_img_sld,'Value',1,'Max',ceil(imgData.t(end)/imgSld_w),'Enable','off')
% set XLim of axes
val_l = imgData.t(1);
val_u = imgData.t(end);
h_ax = [hObjs.ax_img,hObjs.ax_img_sparks,hObjs.ax_prof];
arrayfun(@(x) set(x,'XLim',[val_l val_u]),h_ax)
hObjs.h_bg_sld.SelectedObject = hObjs.h_b1_sld;

% remove listeners
delete(getappdata(mainFig,'propListener'))
delete(getappdata(mainFig,'editScaleListener'))

% remove main figure functions
set(mainFig,'WindowKeyPressFcn',[]);
set(mainFig,'WindowButtonUpFcn',[]);
set(mainFig,'WindowButtonMotionFcn',[]);

set(mainFig,'Pointer','arrow')
drawnow


end

