function reAnalyze(hO,~,mainFig)
% re analyze image, load parameters from old analyis (xls)

% open image
openImg([],[],mainFig)

% get data 
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');

% try to load old analysis (xls file)
% find path to selected data analysis
[~,nImg,~] = fileparts(imgData.fileName);

% look for excel file associated with this image
listing = dir(imgData.filePath);
listing = struct2table(listing);
indXLS = cellfun(@(x) ~isempty( regexp(x,'.xls|.xlsx','once') ), ...
    listing.name);
fileNameXLS = listing.name(indXLS);

indImgXLS = cellfun(@(x) ~isempty( regexp(x,nImg,'once') ), fileNameXLS);
try
%     oldAnalysisInfo = readtable(fullfile(imgData.filePath,fileNameXLS{indImgXLS}),...
%         'ReadVariableNames',false,...
%         'Sheet','info');
%     oldAnalysisInfo = table2cell(oldAnalysisInfo);
    oldAnalysisInfo = readcell( ...
        fullfile(imgData.filePath,fileNameXLS{indImgXLS}), ...
        'Sheet','info');
catch
    disp('did not work, no xls file')
    return
end

% get info from old analysis
[rA,~] = find(cellfun(@(x) ~isempty(strfind(x,'animal:')), ...
    oldAnalysisInfo, 'UniformOutput',1)==1);
[rBlank,~] = find(cellfun(@(x) ~isempty(strfind(x,'blank:')), ...
    oldAnalysisInfo, 'UniformOutput',1)==1);
rPxSzX = find(cellfun(@(x) ~isempty(strfind(x,'px size x (um)')), ...
    oldAnalysisInfo, 'UniformOutput',1)==1);
rPxSzT = find(cellfun(@(x) ~isempty(strfind(x,'px size t (ms)')), ...
    oldAnalysisInfo, 'UniformOutput',1)==1);

[rExpS,~] = find(cellfun(@(x) ~isempty(strfind(x,'exp_notes:')), ...
    oldAnalysisInfo, 'UniformOutput',1)==1);
rExpE = find( cellfun(@(x) all(ismissing(x)), ...
    oldAnalysisInfo(rExpS+1:end,1)), 1, 'first');


animal = oldAnalysisInfo{rA,2};
if isempty(rExpE)
    expInfo = oldAnalysisInfo(rExpS+1:end, 1);
else
    expInfo = oldAnalysisInfo(rExpS+1:rExpS+rExpE-1, 1);
end
blank = oldAnalysisInfo{rBlank,2};

% set up main window
% set up animal
hObjs.popUpMenuAnimal.Value = ...
    find(strcmp(hObjs.popUpMenuAnimal.String,animal));

% set up notes
set(hObjs.h_table_notes,'Data',expInfo)

% set up blank value
hObjs.h_edit_Blank.String = blank;

% set up pizel values in edit fields
try
    hObjs.h_edit_pxSzX.String = num2str(oldAnalysisInfo{rPxSzX,2});
    hO_sim.Tag = 'pxSzX';
    hO_sim.String = num2str(oldAnalysisInfo{rPxSzX,2});
    setPxSzManually(hO_sim, [], mainFig)

    hObjs.h_edit_pxSzT.String = num2str(oldAnalysisInfo{rPxSzT,2});
    hO_sim.Tag = 'pxSzT';
    hO_sim.String = num2str(oldAnalysisInfo{rPxSzT,2});
    setPxSzManually(hO_sim, [], mainFig)
catch 
end

end

