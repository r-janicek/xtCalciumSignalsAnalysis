function tableEditEvents(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% do not let to change numbers only delete
tbl = varargin{1,1};
eventData = varargin{1,2};
mainFig = varargin{1,3};
idx = eventData.Indices(1,1);
selectedROIs = getappdata(mainFig,'selectedROIs');
roiNameToDelete = selectedROIs.roiName(idx);

if any(isnan(eventData.NewData)) || isempty(eventData.NewData)

    % remove row from data table
    data = get(tbl,'Data');
    data(idx,:) = [];
    % delete rectangle and text affiliated with removed object
    delete(selectedROIs.h_rectROI(idx));
    delete(selectedROIs.h_textROI(idx));
    % delete row from selected ROIs
    selectedROIs(idx,:) = [];
 
    % delete event analysis figure, if exists
    % get handles of figures of individial events analysis
    h_figsEventsAnalysis = findall(0,'Type','figure','Tag','imgsOfEventsFromAnalysis');
    figNames = arrayfun(@(x) x.Name,h_figsEventsAnalysis,'UniformOutput',0);
    m_figsEventsAnalysis = strcmp(figNames,roiNameToDelete);
    if any(m_figsEventsAnalysis)
       delete(h_figsEventsAnalysis(m_figsEventsAnalysis)) 
    end
    
else
    % do not let to change string
    data = get(tbl,'Data');
    data(idx,1) = eventData.PreviousData;
      
end

% save changed data
set(tbl,'Data',data);
setappdata(mainFig,'selectedROIs',selectedROIs)

end

