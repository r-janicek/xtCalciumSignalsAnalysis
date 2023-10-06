function getPos(~,~,mainFig)
% get data related to roi of selected profile
% get data
profileAnalysis = getappdata(mainFig,'profileAnalysis');
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');

% check if numbers of detected peaks is the same as number detected events
if numel(profileAnalysis.statEventsSpRec) ~= ...
        numel(profileAnalysis.h_PeaksCirc)   
    errordlg('Not equals numbers of selected events and fits')
    return
end 

n_px = profileAnalysis.numOfPxAvrg;
croppedDataProfile = profileAnalysis.croppedDataProfile;
croppedDataProfileR = profileAnalysis.croppedDataProfileR;
croppedDataWholeAreaRoi = profileAnalysis.croppedDataWholeAreaRoi;
croppedDataWholeAreaRoiRaw = profileAnalysis.croppedDataWholeAreaRoiR;
startOfSpark = profileAnalysis.startOfSpark;        % positions in pixels
endOfSpark = profileAnalysis.endOfSpark;            % positions in pixels
statEventsSpRec = profileAnalysis.statEventsSpRec;  % get properties of selected regions/spark
detectedEventsMask = profileAnalysis.detectedEventsMask;

h_rect_prof = profileAnalysis.h_rect_prof;
h_rect_pos = h_rect_prof.Position;
% h_rect_pos = [floor(h_rect_pos(1)) round(h_rect_pos(2)) h_rect_pos(3) h_rect_pos(4)];
% h_rect_prof.setPosition(h_rect_pos)

ax_img = hObjs.ax_img;
h_PeaksCirc = profileAnalysis.h_PeaksCirc;

pos_ROI_table = h_rect_pos(2) + (h_rect_pos(4)-1)/2;
posOfEvents = profileAnalysis.posOfEvents;
prevFitCoef = profileAnalysis.fitCoefSparkRise;
smooth_span = str2double(get(hObjs.h_smooth_edit,'String'));
bs_crit = str2double(get(hObjs.h_bsDet_edit,'String'));

% plot patch object to show selected profile area
X = [0 h_rect_pos(3) h_rect_pos(3) 0];
Y = [pos_ROI_table+((n_px-1)/2) 
     pos_ROI_table+((n_px-1)/2) 
     pos_ROI_table-((n_px-1)/2) 
     pos_ROI_table-((n_px-1)/2)];
h_patch = patch(X, Y, 'black', 'Parent',ax_img, 'FaceAlpha',0.25);

if isfield(profileAnalysis,'selectedROIs')
    
    T_old = profileAnalysis.selectedROIs;
    % remove result of analysis if any
    m = strcmp(T_old.Properties.VariableNames, 'AnalysisResult');
    m(find(m==1):end) = true;
    T_old(:,m) = [];
    % check if there is a variable whole profile fit
    if any(strcmp(T_old.Properties.VariableNames, 'wholeProfileFit'))
        T_new = table(pos_ROI_table, {croppedDataProfile}, ...
            {croppedDataProfileR}, {croppedDataWholeAreaRoi}, ...
            {croppedDataWholeAreaRoiRaw}, {startOfSpark}, {endOfSpark},...
            {statEventsSpRec}, {detectedEventsMask}, h_patch, ...
            {posOfEvents}, {prevFitCoef}, smooth_span, bs_crit, {nan},...
            'VariableNames',{'positionOfRoi' 'dataProfile' ...
            'dataProfileRaw' 'dataWholeArea' 'dataWholeAreaRaw' ...
            'startPosOfSpark' 'endPosOfSpark' 'regionsProp' ...
            'detectedEventsMask' 'patch' 'eventsPeaks' ...
            'prevFitCoef' 'smooth_span' 'bs_crit' 'wholeProfileFit'});
    else
        T_new = table(pos_ROI_table, {croppedDataProfile}, ...
            {croppedDataProfileR}, {croppedDataWholeAreaRoi}, ...
            {croppedDataWholeAreaRoiRaw}, {startOfSpark}, {endOfSpark},...
            {statEventsSpRec}, {detectedEventsMask}, h_patch, ...
            {posOfEvents}, {prevFitCoef}, smooth_span, bs_crit,...
            'VariableNames',{'positionOfRoi' 'dataProfile' ...
            'dataProfileRaw' 'dataWholeArea' 'dataWholeAreaRaw' ...
            'startPosOfSpark' 'endPosOfSpark' 'regionsProp' ...
            'detectedEventsMask' 'patch' 'eventsPeaks' ...
            'prevFitCoef' 'smooth_span' 'bs_crit'});        
    end
    selectedROIs = [T_old;T_new];
    selectedROIs = sortrows(selectedROIs,{'positionOfRoi'},{'ascend'});  
else
       
    selectedROIs = table(pos_ROI_table, {croppedDataProfile}, ...
        {croppedDataProfileR}, {croppedDataWholeAreaRoi}, ...
        {croppedDataWholeAreaRoiRaw}, {startOfSpark}, {endOfSpark},...
        {statEventsSpRec}, {detectedEventsMask}, h_patch, ...
        {posOfEvents}, {prevFitCoef}, smooth_span, bs_crit,...
        'VariableNames',{'positionOfRoi' 'dataProfile' ...
        'dataProfileRaw' 'dataWholeArea' 'dataWholeAreaRaw' ...
        'startPosOfSpark' 'endPosOfSpark' 'regionsProp' ...
        'detectedEventsMask' 'patch' 'eventsPeaks' ...
        'prevFitCoef' 'smooth_span' 'bs_crit'});
    selectedROIs.Properties.VariableUnits = {'um' '' '' '' '' ...
        'pixels' 'pixels' '' '' '' '[F/F0, ms]' '' 'ms' ''};
end

% set data of table in main fig window
ROIs_table = hObjs.h_table_profs;
tab_data = get(ROIs_table,'Data');
N_events = length(h_PeaksCirc);
if sum(tab_data(:))==0 || isempty(tab_data)     
    set(ROIs_table,'Data',[pos_ROI_table, N_events]);
else
    tab_data = [tab_data;[pos_ROI_table, N_events]];
    tab_data = sortrows(tab_data);
    set(ROIs_table,'Data',tab_data);
    
end
% save data
profileAnalysis.selectedROIs = selectedROIs;
setappdata(mainFig,'profileAnalysis',profileAnalysis);

% set position of patch
uistack(h_patch,'down');

end

