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

% look for pdf/pdfs with this image
listing = dir(imgData.filePath);
listing = struct2table(listing);
indXLS = cellfun(@(x) ~isempty( regexp(x,'.xls|.xlsx','once') ), listing.name);
fileNameXLS = listing.name(indXLS);

indImgXLS = cellfun(@(x) ~isempty( regexp(x,nImg,'once') ), fileNameXLS);
try
%     oldAnalysisInfo = readtable(fullfile(imgData.filePath,fileNameXLS{indImgXLS}),...
%         'ReadVariableNames',false,...
%         'Sheet','info', ...
%         'FileType', 'spreadsheet');
%     oldAnalysisInfo = table2cell(oldAnalysisInfo);
    oldAnalysisInfo = readcell(fullfile(imgData.filePath,fileNameXLS{indImgXLS}), ...
        'Sheet','info');
catch
    disp('did not work, no xls file')
    return
end

% get info from old analysis
[rA,~]=find(cellfun(@(x) ~isempty(strfind(x,'animal:')),oldAnalysisInfo,'UniformOutput',1)==1);
[rBlank,~]=find(cellfun(@(x) ~isempty(strfind(x,'blank:')),oldAnalysisInfo,'UniformOutput',1)==1);
[rExpS,~]=find(cellfun(@(x) ~isempty(strfind(x,'exp_notes:')),oldAnalysisInfo,'UniformOutput',1)==1);
rExpE = find( cellfun(@(x) all(ismissing(x)), oldAnalysisInfo(rExpS+1:end,1)),1,'first');

animal = oldAnalysisInfo{rA,2};
expInfo = oldAnalysisInfo(rExpS+1:rExpS+rExpE-1,1);
blank = oldAnalysisInfo{rBlank,2};

% set up main window
% set up animal
% if isnan(animal)
%     animal = 'wt';
% end
hObjs.popUpMenuAnimal.Value = ...
    find(strcmp(hObjs.popUpMenuAnimal.String,animal));
 
% set up notes
set(hObjs.h_table_notes,'Data',expInfo)

% set up blank value
hObjs.h_edit_Blank.String = blank;


end

