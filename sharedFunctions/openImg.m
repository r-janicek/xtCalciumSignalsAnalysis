function openImg(hO, ~,mainFig)
% open image file using OME-Bioformats

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
if (isfield(getappdata(mainFig),'imgData') && reAnalyzeImgFlag) || ...
        batchProcessing % getappdata(mainFig,'loadFlag')==1
    % re-analyze or batch processing
    filePath = getappdata(mainFig,'FilePath');
    fileName = getappdata(mainFig,'FileName');
    if ~iscell(fileName)
        fileName = {fileName}; 
    end
else
    % load single file
    if isfield(getappdata(mainFig),'imgData')
        % path to previously analyzed image
        imgDataOld = getappdata(mainFig,'imgData');
        pathOld = imgDataOld.filePath;
        % get file path       
        [fileName, filePath] = uigetfile( ...
            {'*.oif;*.oib;*.tif;*.pic;*.msr'}, ...
            'Load image',...
            pathOld, ...
            'Multiselect','off');
    else
        % get file path
        [fileName, filePath] = uigetfile( ...
            {'*.oif;*.oib;*.tif;*.pic;*.msr'}, ...
            'Load image',...
            '',...
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

if isfield(getappdata(mainFig), 'imgData') && ...
        isfield(getappdata(mainFig), 'seriesOfImgs')
    rmappdata(mainFig, 'imgData'); 
    rmappdata(mainFig, 'seriesOfImgs'); 
end
% handles to gui objects
hObjs = getappdata(mainFig, 'hObjs');
% clear axes
if isfield(hObjs, 'ax_img_sparks')
    h_ax = [hObjs.ax_img,hObjs.ax_img_sparks,hObjs.ax_prof];
else
    h_ax = [hObjs.ax_img,hObjs.ax_prof];
end
arrayfun(@(x) cla(x), h_ax)
% get image file extension, used to get pixel size
[~, ~, ext] = fileparts(fileName{1,1});

% load data using OME_bioformats
try
    % load all data about image using OME_bioformats
    data = bfopen(fullfile(filePath,fileName{1,1}));
    % get meta data
    metaData = cat(2, cell(data{1,2}.keySet.toArray), ...
        cell(data{1,2}.values.toArray));
    % get scan mode 
    indx_SM = find(~cellfun(@isempty,regexp(metaData(:,1),'\<ScanMod+e','match'))==1,1,'first');
    try
        scanMode = metaData{indx_SM,2};
    catch
        scanMode = [];
    end
    % create metadata structure
    try
        for i = 1:size(metaData,1)
            ImageDescription_Names{i,1} = genvarname(char(metaData{i,1}));
            ImageDescription_val {i,1} = metaData{i,2};
        end
        [~,ia,~] = unique(ImageDescription_Names);
        metaDataS = cell2struct(cell(ImageDescription_val(ia)), cell(ImageDescription_Names(ia)),1);
    catch 
    end
    % get OME meta data
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
    imgDataXT = data{1,1}(:,1);
    imgDataXT = cellfun(@(x) double(imrotate(x,90)), ...
        imgDataXT, 'UniformOutput',0);

    % assign channels and set up popup menus
    chNames = arrayfun(@(x) sprintf('ch #%d',x), ...
        (1:numel(data{1,1}(:,1))), 'UniformOutput',0);
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
    
    % get pixel sizes
    % time dimension (ms)
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
    % spatial dimension (µm)
    pxSzX = double(omeMetaData.getPixelsPhysicalSizeX(0).value);
    if pxSzX>100000 % looks like error in olympus writing metadata
        pxSzX = pxSzX/1000000;
    end
  
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
        
        numOfROIs = find(cellfun(@(x) ~isempty(x), ...
            strfind(fileROIlines,'ROIBaseFileInformation'))==1);
        
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
        scanline = char(omeMD1( ...
            ~cellfun(@isempty, regexp(omeMD1,'Line'))));
        scanline = char( ...
            regexp(scanline, 'X1=".*" X2=".*" Y1=".*" Y2=".*"', 'match'));
        pos_l = sscanf(scanline, 'X1="%f" X2="%f" Y1="%f" Y2="%f"');
        pos_l = pos_l + 1;
    end
    if ~exist('pos_l','var')
        pos_l = [];
    end

catch
    try
        % if fails with bioformats use matlab builtin function
        data = imread(fullfile(filePath,fileName{1,1}));
        if ~iscell(data)
            data = {data};
        end
        % get metadata
        metaData = imfinfo(fullfile(filePath,fileName{1,1}));
        fieldsMetaData = fieldnames(metaData);
        % save data
        imgDataXT = data;
        % rotate image
        imgDataXT = cellfun(@(x) double(imrotate(x,90)), imgDataXT, ...
            'UniformOutput',0);
        imgDataXTfluo = imgDataXT{1,1};
        imgDataXY = [];
        imgDataXTtrans = [];
        % position of line and pixel size
        pos_l = [];
        pxSzT = 1;
        pxSzX = 1;
        % assign channels and set up popup menus
        hObjs = getappdata(mainFig,'hObjs');
        chNames = arrayfun(@(x) sprintf('ch #%d',x), ...
            (1:numel(data)), 'UniformOutput',0);
        chNames = [chNames,{'none'}];
        set(hObjs.popUpMenuChFluo,'String',chNames,'Value',1)
        set(hObjs.popUpMenuChTrans,'String',chNames,'Value',2)
        
    catch exception
        errordlg('Could not open image file!')
    end
end

% create time vector, starting at 0 (ms)
t = linspace(0, (size(imgDataXTfluo,2)-1)*pxSzT, size(imgDataXTfluo,2));
% filter data
[imgDataXTfluo, ImgFiltersUsed] = imgFiltering( ...
    imgDataXTfluo, pxSzT, pxSzX);
% show image in axes
image(imgDataXTfluo, ...
    'YData',[1 size(imgDataXTfluo,1)], ...
    'XData',[0 max(t)],...
    'CDataMapping','scaled','Parent',hObjs.ax_img);
set(hObjs.ax_img, 'XTick',[], 'FontSize',14, 'Tag','img')
set(get(hObjs.ax_img,'Ylabel'), 'String','x (px)', ...
    'FontWeight','bold')
% time profile
prof_t = mean(imgDataXTfluo, 1);
% show time profile
plot(t,prof_t, 'Parent',hObjs.ax_prof, 'Tag','wholeCellProfile');
set(hObjs.ax_prof, 'FontSize',14, 'Tag','profile')
set(get(hObjs.ax_prof,'Xlabel'), 'String','t (ms)', 'FontWeight','bold')
set(get(hObjs.ax_prof,'Ylabel'), 'String','Fluorescence (F)', ...
    'FontWeight','bold')
xlim(hObjs.ax_prof, [min(t) max(t)])
ylim(hObjs.ax_prof, getAxisLimits(prof_t, 5))
% link axes
try
    linkaxes([hObjs.ax_img, hObjs.ax_img_sparks, hObjs.ax_prof], 'x')
catch
    linkaxes([hObjs.ax_img, hObjs.ax_prof], 'x')
end
% set pixel sizes in edit boxes
hObjs.h_edit_pxSzT.String = num2str(pxSzT);
hObjs.h_edit_pxSzX.String = num2str(pxSzX);

% get axes limits
z_max_img = prctile(imgDataXTfluo(:), 99, 'all');
z_min_img = prctile(imgDataXTfluo(:), 1, 'all');

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
                 'scanLinePos',pos_l,...
                 'microscopeVer',microscopeVer,...
                 'userName',userName);
% save image data in gui data
setappdata(mainFig,'imgData',imgData)
% set up main window
set(hObjs.txt_name_img, ...
    'String',sprintf('%s',fullfile(filePath,char(fileName))),...
    'Tooltip',sprintf('%s',fullfile(filePath,char(fileName))));
% show trans image in small window
image(imgDataXTtrans, 'CDataMapping','scaled', ...
    'Parent',hObjs.h_ax_transCh,...
    'ButtonDownFcn',{@axesTransChFunction});
set(hObjs.h_ax_transCh, 'XTick',[], 'YTick',[], 'PickableParts','all')
set(hObjs.h_ax_transCh, 'ButtonDownFcn',{@axesTransChFunction})

switch getappdata(mainFig,'analysisType')
    case 'spark detection'
        set(hObjs.txt_SpF, ...
            'String',['# sp*100',char(956),'m-1*s-1'])
        set(hObjs.txt_correctedSpF, ...
            'String',{'#';['sp*100',char(956),'m-1*s-1']})
        set(hObjs.h_push_sparks, 'Enable','on')
    
    case 'transients & waves'
        % set normalization panel
        set(hObjs.h_pb_norm,'Enable','on')
        set(hObjs.h_pb_norm,'String','normalize data')
        set(hObjs.h_pb_norm,'Callback',{@norm_data,mainFig})
        set(hObjs.h_pb_norm,'FontWeight','normal')
        set(hObjs.check_doBsFit,'Enable','on')
        % check if baseline fitting is on
        checkboxFcn(hObjs.check_doBsFit, [], mainFig)
        % setup function to select baseline
        set(hObjs.ax_prof, 'buttondownfcn',@mouseSetMaskFcn)
        set(hObjs.h_push_eventImgROI, ...
            'String','<html> <p align="center"> ROI to select image <br> of event <html>')
        set(hObjs.h_push_eventImgROI, 'FontWeight','normal', ...
            'Callback',{@selectEventImage,mainFig})
        % restart table of selected events
        set(hObjs.h_table_eventsImgs,'Data',repmat({'----'},[5,1]))

    otherwise
        set(hObjs.h_table_profs,'Data',zeros(5,2))
        % set normalization panel
        set(hObjs.h_pb_norm,'Enable','on')
        set(hObjs.h_pb_norm,'String','normalize data')
        set(hObjs.h_pb_norm,'Callback',{@norm_data,mainFig})
        set(hObjs.h_pb_norm,'FontWeight','normal')
        set(hObjs.h_push_sparks, 'Enable','off')
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

% set blank panel
set(hObjs.h_push_BlankROI, ...
    'String','set ROI for blank', 'FontWeight','normal')
set(hObjs.h_push_BlankROI, 'Callback',{@setBlankROI,mainFig})
set(hObjs.h_push_BlankROI, 'Enable','on')

set(hObjs.h_edit_Blank, 'String','190') % 'nan' 190 152)
set(hObjs.h_edit_Blank, 'Enable','on')

if isempty(imgDataXY)    
    set(hObjs.h_pb_GetBlank, 'Enable','off')    
else
    set(hObjs.h_pb_GetBlank, 'Enable','on')   
end

% set up pushbutton for cropping
set(hObjs.h_pb_crop, ...
    'String','<html> <p align="center"> set ROI to crop <html>', ...
    'FontWeight','normal', 'Enable','on')
if strcmp(get(hObjs.h_pb_crop,'String'),'crop')  
    set(hObjs.h_pb_crop,'String','set ROI for cropping', ...
        'FontWeight','normal')
end

% set FWHM for spark detection, in pixels, 2 um
try
    set(hObjs.h_edit_fFWHM,'String','10')
    % set(hObjs.h_edit_fFWHM,'String',num2str(ceil(2/pxSzX)))
catch
end
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

if isfield(getappdata(mainFig),'selectedROIs')
    rmappdata(mainFig,'selectedROIs')
end

% restore image slider
imgSld_w = str2double(hObjs.h_edit_imgSld_w.String)*1000; % change from s to ms
set(hObjs.h_img_sld, ...
    'Value',1, ...
    'Max',ceil(imgData.t(end)/imgSld_w), ...
    'Enable','off')
% set XLim of axes
val_l = imgData.t(1);
val_u = imgData.t(end);
arrayfun(@(x) set(x,'XLim',[val_l val_u]), h_ax)
hObjs.h_bg_sld.SelectedObject = hObjs.h_b1_sld;

% remove listeners
delete(getappdata(mainFig,'propListener'))
delete(getappdata(mainFig,'editScaleListener_img'))
delete(getappdata(mainFig,'editScaleListener_sparks'))

% remove main figure functions
set(mainFig,'WindowKeyPressFcn',[]);
set(mainFig,'WindowButtonUpFcn',[]);
set(mainFig,'WindowButtonMotionFcn',[]);

% add listener for interactive scale for image axes
% set and add interactive scale
sc_num = (diff(hObjs.ax_img.YLim)*pxSzX)/2; % in um
hObjs.h_txt_scale_img.String = [sprintf('%0.2f',sc_num),' \mum'];
editScaleListener_img = addlistener(hObjs.ax_img, 'YLim', 'PostSet', ...
    @(varargin)scaleLineChange(varargin,mainFig));
setappdata(mainFig,'editScaleListener_img',editScaleListener_img)

set(mainFig,'Pointer','arrow')
drawnow


end

