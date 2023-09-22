function openImg(~, ~,mainFig)

% get path of image
if isfield(getappdata(mainFig),'reAnalyzeImg')
    if getappdata(mainFig,'reAnalyzeImg')==1
        reAnalyzeImgFlag = true;
    else
        reAnalyzeImgFlag = false;
    end
else    
    reAnalyzeImgFlag = false;    
end

if isfield(getappdata(mainFig),'imgData') && reAnalyzeImgFlag % getappdata(mainFig,'loadFlag')==1
     
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

[~,~,ext] = fileparts(fileName{1,1});

%load all data about image using OME_bioformats
data = bfopen([filePath,fileName{1,1}]);

%meta data
metaData = cat(2,cell(data{1,2}.keySet.toArray),cell(data{1,2}.values.toArray));

indx_SM = find(~cellfun(@isempty,regexp(metaData(:,1),'\<ScanMod+e','match'))==1,1,'first');

try
    scanMode = metaData{indx_SM,2};
catch
    scanMode = [];
end

for i = 1:size(metaData,1)
    
    ImageDescription_Names{i,1} = genvarname(char(metaData{i,1}));
    ImageDescription_val {i,1} = metaData{i,2};
    
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

pxSzX = double(omeMetaData.getPixelsPhysicalSizeX(0).value);
if pxSzX > 10
    pxSzX = pxSzX/1000000;
end

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
    
    pos_l = sscanf(scanline, 'Line X1="%f" X2="%f" Y1="%f" Y2="%f"');
    pos_P = sscanf(pointTP, 'Point X="%f" Y="%f"');
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


ax_prof = hObjs.ax_prof;
plot(t,prof_t,'Parent',ax_prof,'Tag','wholeCellProfile');
set(ax_prof,'FontSize',14)
set(get(ax_prof,'Xlabel'),'String','t (ms)','FontWeight','bold')
set(get(ax_prof,'Ylabel'),'String','Fluorescence (F)','FontWeight','bold')
xlim(ax_prof,[min(t) max(t)])
% mouse baseline selection tool
set(ax_prof,'buttondownfcn',@mouseSetMaskFcn)

try
    ylim(ax_prof,[min(prof_t) max(prof_t)])
catch
    ylim(ax_prof,[min(prof_t) max(prof_t)+1])
end

% link axes
linkaxes([ax_img,ax_prof],'x')

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
else
   
   TPP_delays = [nan;nan];
   posOfTPPinScanLine = [];
   durOfTPP = [];
   s_TPP = [];
   e_TPP = [];
   d_TP_laser = 0;
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
                 'posOfTPPinScanLine',posOfTPPinScanLine);

setappdata(mainFig,'imgData',imgData)

% set up main window
set(hObjs.txt_name_img,'String',sprintf('%s',[filePath,char(fileName)]),...
    'TooltipString',sprintf('%s',[filePath,char(fileName)]));
set(hObjs.h_table_eventsImgs,'Data',repmat({'----'},[5,1]))
% show trans image in small window
image(imgDataXTtrans,'CDataMapping','scaled','Parent',hObjs.h_ax_transCh,...
    'ButtonDownFcn',{@axesTransChFunction});
set(hObjs.h_ax_transCh,'XTick',[],'YTick',[],'PickableParts','all')
set(hObjs.h_ax_transCh,'ButtonDownFcn',{@axesTransChFunction})

% set normalization panel
set(hObjs.h_pb_norm,'Enable','on')
set(hObjs.h_pb_norm,'String','normalize data')
set(hObjs.h_pb_norm,'Callback',{@norm_data,mainFig})
set(hObjs.h_pb_norm,'FontWeight','normal')

set(hObjs.check_doBsFit,'Enable','on','Value',0)

% set blank panel
set(hObjs.h_push_BlankROI,'String','set ROI for blank','FontWeight','normal')
set(hObjs.h_push_BlankROI,'Callback', {@setBlankROI,mainFig})
set(hObjs.h_push_BlankROI,'Enable','on')

set(hObjs.h_edit_Blank,'String','nan')
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

% keyboard
% hObjs.h_table_notes.Data = [hObjs.h_table_notes.Data; {''}];

% % remove old analysis
% if isfield(getappdata(mainFig),'lastPressedPusbutton')
%     rmappdata(mainFig,'lastPressedPusbutton');   
% end
% 
% if isfield(getappdata(mainFig),'profileAnalysis')
%     rmappdata(mainFig,'profileAnalysis');   
% end

% remove old selected ROIs
if isfield(getappdata(mainFig),'selectedROIs')
    rmappdata(mainFig,'selectedROIs');   
end

% delete any hidden windows from analysis
% get handles of figures of individial events analysis
delete(findall(0,'Type','figure','Tag','imgsOfEventsFromAnalysis'));

% restore image slider
imgSld_w = str2double(hObjs.h_edit_imgSld_w.String)*1000; % change from s to ms
set(hObjs.h_img_sld,'Value',1,'Max',ceil(imgData.t(end)/imgSld_w),'Enable','off')
% set XLim of axes
val_l = imgData.t(1);
val_u = imgData.t(end);
h_ax = [hObjs.ax_img,hObjs.ax_prof];
arrayfun(@(x) set(x,'XLim',[val_l val_u]),h_ax)
hObjs.h_bg_sld.SelectedObject = hObjs.h_b1_sld;

% remove listeners
delete(getappdata(mainFig,'propListener'))
delete(getappdata(mainFig,'editScaleListener'))

% add listener
% set and add interactive scale
sc_num = (diff(ax_img.YLim)*pxSzX)/2; % in um
hObjs.h_txt_scale.String = [sprintf('%0.2f',sc_num),' \mum'];
editScaleListener = addlistener(ax_img,'YLim','PostSet',@(varargin)scaleLineChange(varargin,mainFig));
setappdata(mainFig,'editScaleListener',editScaleListener)

% reset buttons for ROIs and table
delete(findobj(ax_img,'Tag','imrect'))
set(hObjs.h_push_eventImgROI,'String','<html> <p align="center"> ROI to select image <br> of event <html>')
set(hObjs.h_push_eventImgROI,'Callback',{@selectEventImage,mainFig})
set(hObjs.h_push_eventImgROI,'FontWeight','normal')

set(hObjs.h_table_eventsImgs,'Data',repmat({'----'},[5,1]))

% remove main figure functions
set(mainFig,'WindowKeyPressFcn',[]);
set(mainFig,'WindowButtonUpFcn',[]);
set(mainFig,'WindowButtonMotionFcn',[]);

set(mainFig,'Pointer','arrow')
drawnow

end

