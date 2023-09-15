function tableEdit(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% do not let to change numbers only delete

tbl = varargin{1,1};
eventData = varargin{1,2};
h_mainFig = varargin{1,3};
idx = eventData.Indices(1,1);
profileAnalysis = getappdata(h_mainFig,'profileAnalysis');

selectedROIs = profileAnalysis.selectedROIs;

if isnan(eventData.NewData)

    % remove data
    data = get(tbl,'Data');
    data(idx,:) = [];
    p_h = selectedROIs.patch(idx);
    
    % delete corresponding patch object
    selectedROIs(idx,:) = [];
    delete(p_h);
          
else
    % do not let to change numbers
    data = get(tbl,'Data');
    data(idx,1) = eventData.PreviousData;
      
end

profileAnalysis.selectedROIs = selectedROIs;

% save changed data
set(tbl,'Data',data);
setappdata(h_mainFig,'profileAnalysis',profileAnalysis)

end

