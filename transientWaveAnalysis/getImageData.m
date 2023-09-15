function getImageData(~,~,mainFig,h_rectROI)

% get data
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');
ax_img = hObjs.ax_img;
hROIsTable = hObjs.h_table_eventsImgs;
eventType = hObjs.h_bg_eventSelection.SelectedObject.String;

% create name of ROI

try
    numIdxW = cellfun(@(x) str2double(x), cellfun(@(x) regexp(x,'\d*','match'), hROIsTable.Data) );
    numIdxW(cellfun(@(x) isempty(regexp(x,'wave','match')), hROIsTable.Data)) = [];
catch
    numIdxW = 0;
end
try
    numIdxTt = cellfun(@(x) str2double(x), cellfun(@(x) regexp(x,'\d*','match'), hROIsTable.Data) );
    numIdxTt(cellfun(@(x) isempty(regexp(x,'trigg. transient','match')), hROIsTable.Data)) = [];
catch
    numIdxTt = 0;
end
try
    numIdxTs = cellfun(@(x) str2double(x), cellfun(@(x) regexp(x,'\d*','match'), hROIsTable.Data) );
    numIdxTs(cellfun(@(x) isempty(regexp(x,'spon. transient','match')), hROIsTable.Data)) = [];
catch
    numIdxTs = 0;
end
try
    numIdxC = cellfun(@(x) str2double(x), cellfun(@(x) regexp(x,'\d*','match'), hROIsTable.Data) );
    numIdxC(cellfun(@(x) isempty(regexp(x,'caff','match')), hROIsTable.Data)) = [];
catch
    numIdxC = 0;
end
posMissingNumW = find(diff(numIdxW)>1);
posMissingNumTt = find(diff(numIdxTt)>1);
posMissingNumTs = find(diff(numIdxTs)>1);
posMissingNumC = find(diff(numIdxC)>1);

if ~isempty(posMissingNumW)
    indW = numIdxW(posMissingNumW)+1;
else
    if ~isempty(numIdxW)
        indW = max(numIdxW)+1;
    else
        indW = 1;
    end
end

if ~isempty(posMissingNumTt)
    indTt = numIdxTt(posMissingNumTt)+1;
else
    if ~isempty(numIdxTt)
        indTt = max(numIdxTt)+1;
    else
        indTt = 1;
    end
end

if ~isempty(posMissingNumTs)
    indTs = numIdxTs(posMissingNumTs)+1;
else
    if ~isempty(numIdxTs)
        indTs = max(numIdxTs)+1;
    else
        indTs = 1;
    end
end

if ~isempty(posMissingNumC)
    indC = numIdxC(posMissingNumC)+1;
else
    if ~isempty(numIdxC)
        indC = max(numIdxC)+1;
    else
        indC = 1;
    end
end

switch eventType
    case 'wave'
        roiName = sprintf('wave #%d',indW);
        
    case 'transient (triggered)'
        roiName = sprintf('trigg. transient #%d',indTt);
        
    case 'transient (spon.)'
        roiName = sprintf('spon. transient #%d',indTs);
        
    case 'caffeine'
        roiName = sprintf('caffeine #%d',indC);
end
    
% get mask of ROI and cropp data
m_ROI = createMask(h_rectROI);
[r,c] = find(m_ROI);

imgROI_F = imgData.imgDataXTfluoF(min(r):max(r),min(c):max(c)); % filtred
imgROI_FN = imgData.imgDataXTfluoFN(min(r):max(r),min(c):max(c)); % filtred and normalized

t = (0:size(imgROI_FN,2)-1).*imgData.pxSzT;
x = (0:size(imgROI_FN,1)-1).*imgData.pxSzX;

posROI = {getPosition(h_rectROI)};
dataROI = struct('imgF',imgROI_F,...
                 'imgFN',imgROI_FN,...
                 't',t,...
                 'x',x);
         
% color for rectangle
% 32 shades
c = linspace(0,0.9,32);
c = c(:);
switch eventType
    case 'wave'
        recColor = [ones(size(c)) c c];  
 
    case 'transient (triggered)'
        recColor = [c c c];  
   
    case 'transient (spon.)'
        recColor = [ones(size(c)) ones(size(c)) c];
        
    case 'caffeine'
        recColor = [c ones(size(c)) c];
end

ind_color = randi([1 32],1,1);
      
% plot selected ROI, also with text
xp = [posROI{1}(1) posROI{1}(1)+posROI{1}(3) posROI{1}(1)+posROI{1}(3) posROI{1}(1)];
yp = [posROI{1}(2) posROI{1}(2) posROI{1}(2)+posROI{1}(4) posROI{1}(2)+posROI{1}(4)];
rectOfROI = patch(ax_img,'XData',xp,'YData',yp,'FaceColor',recColor(ind_color,:),...
    'FaceAlpha',0.15,'EdgeColor',recColor(ind_color,:),'LineWidth',3);
%rectangle(ax_img,'Position',cell2mat(posROI),'EdgeColor',recColor(ind_color,:),'LineWidth',4);

h_text = text(posROI{1}(1),posROI{1}(2),roiName,'Parent',ax_img,...
        'FontSize',16,'FontWeight','bold','VerticalAlignment','top',...
        'Color',recColor(ind_color,:));       
             
% save ROIs data
if isfield(getappdata(mainFig),'selectedROIs')
  
    T_old = getappdata(mainFig,'selectedROIs');        
    T_new = table({roiName},posROI,dataROI,rectOfROI,h_text,{[]},...
                 'VariableNames',{'roiName' 'positionOfRoi' 'dataROIs' ...
                 'h_rectROI' 'h_textROI' 'analysis'});
    
    selectedROIs = [T_old;T_new];
      
else
       
    selectedROIs = table({roiName},posROI,dataROI,rectOfROI,h_text,{[]},...
                 'VariableNames',{'roiName' 'positionOfRoi' 'dataROIs' ...
                 'h_rectROI' 'h_textROI' 'analysis'});

    selectedROIs.Properties.VariableUnits = {'' '' '' '' '' ''};
    
end

% sort selectedROIs
numIdx = cellfun(@(x) str2double(x), cellfun(@(x) regexp(x,'\d*','match'), selectedROIs.roiName) );
[~,I] = sort(numIdx);
selectedROIs = selectedROIs(I,:);

%selectedROIs = sortrows(selectedROIs,{'roiName'},{'ascend'});

% show in table
set(hROIsTable,'Data',selectedROIs.roiName);

% delete ROI, reset pushbuttons
delete(h_rectROI)
set(hObjs.h_push_eventImgROI,'String','<html> <p align="center"> ROI to select image <br> of event <html>')
set(hObjs.h_push_eventImgROI,'Callback',{@selectEventImage,mainFig})
set(hObjs.h_push_eventImgROI,'FontWeight','normal')

% enable change of type of event
hObjs.h_b1_eventSelection.Enable = 'on';
hObjs.h_b2_eventSelection.Enable = 'on';
hObjs.h_b3_eventSelection.Enable = 'on';
hObjs.h_b4_eventSelection.Enable = 'on';

% save data
setappdata(mainFig,'selectedROIs',selectedROIs)

end

