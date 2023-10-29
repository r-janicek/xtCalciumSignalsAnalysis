function plotROIprofile(profROI, ~, mainFig)

set(mainFig,'Pointer','watch')
drawnow
% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
profileAnalysis = getappdata(mainFig,'profileAnalysis');
analysisType = getappdata(mainFig,'analysisType');
% get image data
img = imgData.imgDataXTfluoFN;
% px size
pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;
t = imgData.t;

ax_img_sparks = hObjs.ax_img_sparks;
crit = str2double(get(hObjs.h_edit_tresh_prof,'String'));
smooth_span = str2double(get(hObjs.h_smooth_edit,'String'));
bs_crit = str2double(get(hObjs.h_bsDet_edit,'String'));
ax_prof = hObjs.ax_prof;

% create ROI mask
profROI_m = createMask(profROI);
r_t = find(any(profROI_m, 1));
r_x = find(any(profROI_m, 2));
% crop from data for selected analysis
cropRoi = img(r_x,:);
% crop from raw data
cropRoiR = imgData.imgDataXTfluoR(r_x,:); 

y_px = linspace(1,size(cropRoi,1),size(cropRoi,1));
y_um = linspace(0,(size(cropRoi,1)-1)*pxSzX,size(cropRoi,1));
c = round(y_px(end)/2);

% calculate average from area +- h_d/2 um around centre (line)
h_d = str2double(get(hObjs.h_edit_averageWidth,'String'));

n_px = ceil(h_d/pxSzX);
if mod(n_px,2)==0    
    n_px = n_px - 1;   
end

% image of cropped ROI 
% image(cropRoi,'YData',[y_um(1) y_um(end)],'XData',[t(1) t(end)],'CDataMapping','scaled','Parent',ax_img_sparks_2);
% set(get(ax_img_sparks_2,'Ylabel'),'String','x (\mum)','FontWeight','bold')
% set(ax_img_sparks_2,'XTick',[],'FontSize',14,'YDir','reverse','YAxisLocation','right')
%keyboard
if strcmp(hObjs.h_bg_sld.SelectedObject.String,'#s part')
    % get old axis limits
    oldXlim = get(ax_img_sparks,'XLim');
end

image(cropRoi, 'YData',[y_px(1) y_px(end)], ...
    'XData',[t(1) t(end)],...
    'CDataMapping','scaled', ...
    'Parent',ax_img_sparks);
set(get(ax_img_sparks,'Ylabel'), 'String','x (pixels)', ...
    'FontWeight','bold', 'Tag','img_sparks')
set(ax_img_sparks,'XTick',[], 'FontSize',14, 'YDir','reverse', ...
    'YAxisLocation','left')

% set and add interactive scale
sc_num = (diff(ax_img_sparks.YLim)*pxSzX)/2; % in um
hObjs.h_txt_scale_sparks.String = [sprintf('%0.2f',sc_num),' \mum'];
editScaleListener = addlistener(ax_img_sparks, 'YLim', 'PostSet', ...
    @(varargin)scaleLineChange(varargin,mainFig));
setappdata(mainFig,'editScaleListener_sparks',editScaleListener)

% show borders of area taken for average profile
line([t(1) t(end)],[c c], 'Parent',ax_img_sparks, ...
    'LineWidth',3, 'Color','k', 'LineStyle','-')
line([t(1) t(end)],[c-((n_px-1)/2) c-((n_px-1)/2)], ...
    'Parent',ax_img_sparks, 'LineWidth',2, 'Color','k', 'LineStyle',':')
line([t(1) t(end)],[c+((n_px-1)/2) c+((n_px-1)/2)], ...
    'Parent',ax_img_sparks, 'LineWidth',2, 'Color','k', 'LineStyle',':')

% standard analysis
% save cropped image data
profileAnalysis.croppedDataWholeAreaRoi = cropRoi; % filtred and possibly normalized data
profileAnalysis.croppedDataWholeAreaRoiR = cropRoiR; % raw data
profileAnalysis.numOfPxAvrg = n_px;
setappdata(mainFig,'profileAnalysis',profileAnalysis)
%find individuals events
try
    statEvents = findEvents(mainFig, cropRoiR);
catch
    statEvents = [];
end
% remove events which bounding boxes are not crossing center lines
try
    lines_px_r = (c-((n_px-1)/2):c+((n_px-1)/2));
    events_px_r = arrayfun(@(x) x.SubarrayIdx{1,1}, statEvents, ...
        'UniformOutput',0);
    m_center = cellfun(@(x) any(ismember(lines_px_r,x)), events_px_r, ...
        'UniformOutput',1);
    statEvents(~m_center) = [];
catch
    statEvents = [];
end

% find centre of the event
pos_centre_x = arrayfun(@(x) x.WeightedCentroid(1,1), statEvents, ...
    'UniformOutput',1);
pos_centre_y = round(arrayfun(@(x) x.WeightedCentroid(1,2), statEvents, ...
    'UniformOutput',1));
col = zeros(length(statEvents), 1);
for i = 1:length(statEvents)
    event = cropRoi(statEvents(i).PixelIdxList);
    [~,p_m] = max(event);
    px = statEvents(i).PixelIdxList;

    col(i,1) = ceil(px(p_m)/size(cropRoi,1));
    % row(i,1) = rem(px(p_m),size(cropRoi,1));
end

try
    pos_max = num2cell([pos_centre_y, col],2);
catch
    pos_max = [];
end

% show centre and boundary rectangle of detected events
prof_t_evnts_m = false(size(cropRoi,2),1);
for i=1:length(statEvents)
    % bounding rectangle
    pos = statEvents(i).BoundingBox;
    eventRec = rectangle( ...
        'Position',[pos(1)*pxSzT pos(2) pos(3)*pxSzT pos(4)],...
        'Parent',ax_img_sparks, 'EdgeColor','r', 'LineStyle',':',...
        'LineWidth',2, 'Tag',num2str(i));
    % centre line
    centre = pos_max{i};
    centreLine = line( ...
        [centre(2)*pxSzT-pxSzT centre(2)*pxSzT-pxSzT], ...
        [centre(1) centre(1)],...
        'Parent',ax_img_sparks, 'Color','r', 'LineStyle','none',...
        'Marker','+', 'MarkerSize',20, 'LineWidth',3, ...
        'PickableParts','all',...
        'ButtonDownFcn',{@profileROIButtonDownFcn, mainFig});
    [statEvents(i).centreLine] = centreLine;
    [statEvents(i).eventRec] = eventRec;
    % create events mask for profile
    prof_t_evnts_m(statEvents(i).SubarrayIdx{2}, 1) = true;
end

% take area for averaging and calculate time profile
cropProf = cropRoi(c-((n_px-1)/2):c+((n_px-1)/2),:);
prof_t(:,1) = mean(cropProf, 1);
cropProfR = cropRoiR(c-((n_px-1)/2):c+((n_px-1)/2),:);
profileAnalysis.croppedDataProfile = cropProf;
profileAnalysis.croppedDataProfileR = cropProfR;

% smooth profile with loess with defined duration in ms
n_pts = round(smooth_span/pxSzT);
prof_t_s = smooth(prof_t,n_pts/length(prof_t),'loess');

% find peaks
noise_std  = std(prof_t(prof_t < (mean(prof_t)+2*std(prof_t))));
mph = mean(prof_t)+crit*noise_std; % minimum peak height
mpw = 10; % minimum peak width in ms
mpd = mpw/2;% minimum peak distance
mpp = 2*noise_std; % minimum peak prominence

if max(prof_t) > mph
    [pks,locs] = findpeaks(prof_t,t,'MinPeakHeight',mph,...
        'MinPeakWidth',mpw,...
        'MinPeakDistance',mpd,...
        'MinPeakProminence',mpp);
end

% plot profile
pp = plot(t,prof_t, t,prof_t_s, 'Parent',ax_prof);
set(pp(1),'Color','k')
set(pp(2),'Color','b', 'LineStyle','-')

if exist('oldXlim', 'var')
    xlim(ax_prof, [oldXlim(1) oldXlim(2)])
else
    xlim(ax_prof, [min(t) max(t)])
end
ylim(ax_prof, getAxisLimits(prof_t,5))
set(ax_prof,'FontSize',14)
% set callback function
set([ax_prof; pp], 'buttondownfcn', {@profileROIButtonDownFcn,mainFig})

% scale size of circles
scSz = get(0,'screensize');
cSz = (min(scSz(3),scSz(4))/1440)*60;
% plot detected peaks in profile
if exist('pks','var')
    if ~isempty(pks)
        h_circ = zeros(length(pks),1);
        for i=1:length(pks)
            % plot detected peaks
            h_circ(i,1) = line(locs(i), pks(i), ...
                'Parent',ax_prof, 'Color','r',...
                'LineStyle','none', 'Marker','.', ...
                'MarkerSize',cSz, 'PickableParts','all',...
                'ButtonDownFcn',{@profileROIButtonDownFcn,mainFig});
        end
        % fit sparks and plot results
        [h_line, detectedEventsMask, coef, ~, startOfSpark, endOfSpark] = ...
            fitSparkRise(pxSzT, t, prof_t, pks, locs, ax_prof, ...
            [], 1e-3, 400, smooth_span, bs_crit, [], [], prof_t_evnts_m);
    end
end

% set correct axis labels
if isfield(imgData,'norm_flag') && ...
        imgData.norm_flag==1 && ...
        strcmp(analysisType,'spark recovery ryanodine')
    set(get(ax_prof,'Xlabel'), 'String','t (ms)', 'FontWeight','bold')
    set(get(ax_prof,'Ylabel'), ...
        'String',['fluorescence ' '(',char(916),'F/F0)'], ...
        'FontWeight','bold')
else
    set(get(ax_prof,'Xlabel'), 'String','t (ms)', 'FontWeight','bold')
    set(get(ax_prof,'Ylabel'), 'String','fluorescence (F)', ...
        'FontWeight','bold')
end

% save data
if exist('h_circ','var') && ~isempty(h_circ)

    profileAnalysis.h_PeaksCirc = h_circ;
    profileAnalysis.posOfEvents = num2cell([pks,locs']);
    profileAnalysis.statEventsSpRec = statEvents;
    profileAnalysis.fitCoefSparkRise = coef;
    profileAnalysis.startOfSpark = startOfSpark;
    profileAnalysis.endOfSpark = endOfSpark;
    profileAnalysis.h_FitLine = h_line;
    profileAnalysis.detectedEventsMask = detectedEventsMask;
    profileAnalysis.peaksCircleSz = cSz;

else
    profileAnalysis.h_PeaksCirc = [];
    profileAnalysis.posOfEvents = num2cell([[],[]]);
    profileAnalysis.statEventsSpRec = statEvents;
    profileAnalysis.fitCoefSparkRise = [];
    profileAnalysis.startOfSpark = [];
    profileAnalysis.endOfSpark = [];
    profileAnalysis.h_FitLine = [];
    profileAnalysis.detectedEventsMask = [];
    profileAnalysis.peaksCircleSz = cSz;
end

setappdata(mainFig,'profileAnalysis',profileAnalysis)



set(mainFig,'Pointer','arrow')
drawnow

end




