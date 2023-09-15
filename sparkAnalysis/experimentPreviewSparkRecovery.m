function experimentPreviewSparkRecovery(~,~,mainFig)

% set up pointer
set(mainFig,'Pointer','watch')
drawnow

% get data 
hObjs = getappdata(mainFig,'hObjs');

if isfield(getappdata(mainFig),'seriesImgData')
    seriesImgData = getappdata(mainFig,'seriesImgData');
    fileNames = seriesImgData.nameImg;      % cell array
    filePath = [seriesImgData.pathCell{1},'/'];     % string
    m_ctrl = seriesImgData.t_s == 0;
    drugN = unique(seriesImgData.drug(~m_ctrl));
    drugC = unique(seriesImgData.c_drug_uM(~m_ctrl));
    t_drugStart = seriesImgData.t_s(find(~m_ctrl,1,'first'));
    
else
    try
        data_setUpFig = getappdata(findobj('Name','set up images'));
        oldPath = data_setUpFig.imgData(1).filePath;
    catch
        try
            data_mainFig = getappdata(mainFig);
            oldPath = data_mainFig.imgData.filePath;
        catch
            oldPath = '/Users/radoslavjanicek/Desktop/analyzedData';
        end
    end
    
    % load time series of images
    [fileNames, filePath] = uigetfile({'*.oif;*.oib;*.tif;*.pic'}, 'Load image',...
        oldPath, 'Multiselect','on');
    
    if ~iscell(fileNames)
        fileNames = {fileNames};
    end
    
end

%% get image data
for i = 1:numel(fileNames)
    imgData(i) = loadImgDataPrev(filePath,fileNames{i});
end

% calculate differences between starts of aquisition of images
dt_imgs = seconds(diff([imgData(:).imgAcquisitionTime])); % in s
t_imgs = [0,cumsum(dt_imgs)];

if exist('t_drugStart','var')
    t0 = t_drugStart; % start of recording of first image with treatment, in s
    % baseline img based on mask
    
    t_imgs_adj = (t_imgs - t_imgs(find(~m_ctrl,1,'first'))) + t0; % time axes adjusted to start of drug aplication
    mBs = m_ctrl;
    
    % prepare data for table
    d = [ fileNames(:), num2cell(mBs), num2cell(t_imgs_adj'),...
        [repmat({'ctrl'},[numel(fileNames(mBs)),1]); repmat(drugN,[numel(fileNames(~mBs)),1])],...
        num2cell( [zeros(numel(fileNames(mBs)),1);ones(numel(fileNames(~mBs)),1).*drugC] ) ];
    
else
    t0 = 30; % start of recording of first image with treatment, in s
    % consider 1. and 2. images as baseline
    if numel(fileNames)<3
        t_imgs_adj = (t_imgs - t_imgs(2)) + t0; % time axes adjusted to start of drug aplication
        mBs = [true;false(numel(fileNames)-1,1)];
    else
        t_imgs_adj = (t_imgs - t_imgs(3)) + t0; % time axes adjusted to start of drug aplication
        mBs = [true;true;false(numel(fileNames)-2,1)];
    end
    
    % prepare data for table
    drugN = {'PDP3'};
    drugC = 1;
    
    d = [ fileNames(:), num2cell(mBs), num2cell(t_imgs_adj'),...
        [repmat({'ctrl'},[numel(fileNames(mBs)),1]); repmat(drugN,[numel(fileNames(~mBs)),1])],...
        num2cell( [zeros(numel(fileNames(mBs)),1);ones(numel(fileNames(~mBs)),1).*drugC] ) ];
    
    
    
end
%% figure
% close previous opened
delete(findobj('Name','set up images'))

% create new
setUpFig = figure('Name','set up images','units','normalized','outerposition',[0.1 0.5 0.5 0.5]);

% show path to image series
uicontrol('Style','text',...
    'Parent',setUpFig,'Units','normalized','Position', [0.05 0.85 0.9 0.1],...
    'FontUnits','normalized','FontSize',0.3,'FontWeight','bold',...
    'HorizontalAlignment','left','String',sprintf('%s',filePath));

% table
h_tbl = uitable('Data',d,'ColumnEditable',[false true false false false],...
    'ColumnName',{'<HTML> <p align="center"> <font size="5"> <b> _____image name_____ </b> </HTML>',...
                  '<HTML> <p align="center"> <font size="4"> <b> baseline <br/> image </b> </HTML>',...
                  '<HTML> <p align="center"> <font size="4"> <b> start of img REC <br/> after drug addition(s) </b> </HTML>',...
                  '<HTML> <p align="center"> <font size="4"> <b> drug name </b> </HTML>',...
                  '<HTML> <p align="center"> <font size="4"> <b> drug conc (uM) </b> </HTML>'},...
    'ColumnWidth','auto','Parent',setUpFig,...
    'Units','normalized','Position',[0.05 0.15 0.6 0.7],...
    'FontUnits','normalized','FontSize',0.05,...
    'CellEditCallback',{@tableEditExpPrev,setUpFig});

% some interaction
uicontrol('Style','text',...
    'Parent',setUpFig,'Units','normalized','Position', [0.7 0.7 0.15 0.1],...
    'FontUnits','normalized','FontSize',0.225,'FontWeight','bold',...
    'HorizontalAlignment','left','String',{'start of recording of first';'treated image (s)'});
h_t_drug = uicontrol('Style','edit','String',num2str(t0),...
    'Parent',setUpFig,'Units','normalized','Position', [0.85 0.725 0.05 0.075],...
    'FontUnits','normalized','FontSize',0.5,...
    'Callback',{@tableEditExpPrev,setUpFig});

uicontrol('Style','text',...
    'Parent',setUpFig,'Units','normalized','Position', [0.7 0.575 0.125 0.1],...
    'FontUnits','normalized','FontSize',0.225,'FontWeight','bold',...
    'HorizontalAlignment','left','String','drug name:');
h_N_drug = uicontrol('Style','edit','String',drugN{1},...
    'Parent',setUpFig,'Units','normalized','Position', [0.85 0.6 0.1 0.075],...
    'FontUnits','normalized','FontSize',0.5,...
    'Callback',{@tableEditExpPrev,setUpFig});

uicontrol('Style','text',...
    'Parent',setUpFig,'Units','normalized','Position', [0.7 0.45 0.15 0.1],...
    'FontUnits','normalized','FontSize',0.225,'FontWeight','bold',...
    'HorizontalAlignment','left','String',{'c of drug (uM)'});
h_c_drug = uicontrol('Style','edit','String',num2str(drugC),...
    'Parent',setUpFig,'Units','normalized','Position', [0.85 0.475 0.05 0.075],...
    'FontUnits','normalized','FontSize',0.5,...
    'Callback',{@tableEditExpPrev,setUpFig});

% show whole experiment
h_pb_delay = uicontrol('Style', 'pushbutton',...
    'String','<html> <p align="center"> show time <br/> series of images <html>',...
    'FontUnits','normalized','FontSize',0.2,'FontWeight','bold',...
    'Parent',setUpFig,'Units','normalized','Position', [0.725 0.275 0.15 0.15],...
    'Callback',{@showWholeImageTimeSeries,setUpFig,mainFig},'Enable','on');

% save data
setappdata(setUpFig,'imgData',imgData)
setappdata(setUpFig,'h_tbl',h_tbl)
setappdata(setUpFig,'h_t_drug',h_t_drug)
setappdata(setUpFig,'h_N_drug',h_N_drug)
setappdata(setUpFig,'h_c_drug',h_c_drug)

% show immediately, when calling it from analysis
if exist('t_drugStart','var')
    showWholeImageTimeSeries([],[],setUpFig,mainFig)
end

% set up pointer
set(mainFig,'Pointer','arrow')
drawnow

%%

% %%%%%%%%%%%%%%
   function imgData = loadImgDataPrev(filePath,fileName)
      
        [~,~,ext] = fileparts(fileName);
        
        %load all data about image using OME_bioformats
        data = bfopen([filePath,fileName]);
        
        %meta data
        metaData = cat(2,cell(data{1,2}.keySet.toArray),cell(data{1,2}.values.toArray));
        
        indx_SM = find(~cellfun(@isempty,regexp(metaData(:,1),'\<ScanMod+e','match'))==1,1,'first');
        
        try
            scanMode = metaData{indx_SM,2};
        catch
            scanMode = [];
        end
        
        for k = 1:size(metaData,1)
            
            ImageDescription_Names{k,1} = genvarname(char(metaData{k,1}));
            ImageDescription_val {k,1} = metaData{k,2};
            
        end
        
        [~,ia,~] = unique(ImageDescription_Names);
        metaDataS = cell2struct(cell(ImageDescription_val(ia)), cell(ImageDescription_Names(ia)),1);
        
        %OME meta data
        omeMetaData = data{1,4};
        omeMetaDataXML = char(omeMetaData.dumpXML());
        
        omeMD = textscan(omeMetaDataXML,'%s','Delimiter',{'<','>'});
        omeMD1 = omeMD{1,1}(~cellfun(@isempty,regexp(omeMD{1,1},'=')));
        
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
                
        end
        
        % get time of image acquisition
        try
            dateStr = metaDataS.Global0x5BAcquisitionParametersCommon0x5DImageCaputreDate;
        catch
            keyboard
        end
        
        % remove ' '
        dateStr = regexp(dateStr,'''','split');
        dateStr(cellfun(@(x) isempty(x), dateStr)) = [];
        
        % convert string to date
        imgAcquisitionTime = datetime(dateStr{1},'Format','yyyy-MM-dd HH:mm:ss');
        
        % pixel size in x
        pxSzX = double(omeMetaData.getPixelsPhysicalSizeX(0).value);
         
        t = linspace(0,(size(imgDataXTfluo,2)-1)*pxSzT,size(imgDataXTfluo,2)); %ms
        
        % try to load positions of ROIs from image metadata
        if strcmp(ext,'.oif')
            
            fileList = dir([filePath,fileName{1,1},'.files']);
            fNnames = {fileList.name};
            roiFile = fNnames{cellfun(@(x) ~isempty(x),strfind(fNnames,'.roi'))==1};
            roiFile = fullfile([filePath,fileName{1,1},'.files'],roiFile);
            
            fileID = fopen(roiFile);
            fileROIlines = textscan(fileID,'%s','Delimiter','\n');
            fileROIlines = regexprep(fileROIlines{1,1},'[^a-zA-Z_-.,=0-9]','');
            fclose(fileID);
            
            numOfROIs = find(cellfun(@(x) ~isempty(x),strfind(fileROIlines,'ROIBaseFileInformation'))==1);
            
            for j=1:length(numOfROIs)
                
                try
                    ROI = fileROIlines(numOfROIs(j):numOfROIs(j+1)-1);
                catch
                    ROI = fileROIlines(numOfROIs(j):end);
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
            'imgAcquisitionTime',imgAcquisitionTime,...
            't',t,...
            'imgFiltersUsed',ImgFiltersUsed,...
            'z_max_img',z_max_img,...
            'z_min_img',z_min_img,...
            'scanLinePos',pos_l);
        
        
    end


   function tableEditExpPrev(hO,E,setUpFig)
       
       % get data
       nDrug = get(getappdata(setUpFig,'h_N_drug'),'String');
       cDrug = str2double(get(getappdata(setUpFig,'h_c_drug'),'String'));
       % start of recording of first image with treatment, in s
       t0 = str2double(get(getappdata(setUpFig,'h_t_drug'),'String'));
       imgData = getappdata(setUpFig,'imgData');
       h_tbl = getappdata(setUpFig,'h_tbl');
       d_tbl = get(h_tbl,'Data');
       mBs = [d_tbl{:,2}];
       % index of first image with treatment
       ind = find(mBs==false,1,'first');
       
       % calculate differences between starts of aquisition of images
       dt_imgs = seconds(diff([imgData(:).imgAcquisitionTime])); % in s
       t_imgs = [0,cumsum(dt_imgs)];
       
       % time axes adjusted to start of drug aplication
       t_imgs_adj = (t_imgs - t_imgs(ind)) + t0; 
       
       % save data to table
       d_tbl = [ fileNames(:), num2cell(mBs)', num2cell(t_imgs_adj'),...
                [repmat({'ctrl'},[numel(fileNames(mBs)),1]); repmat({nDrug},[numel(fileNames(~mBs)),1])],...
                num2cell([zeros(numel(fileNames(mBs)),1);ones(numel(fileNames(~mBs)),1).*cDrug]) ];
       
       h_tbl.Data = d_tbl;
       
    end


end

