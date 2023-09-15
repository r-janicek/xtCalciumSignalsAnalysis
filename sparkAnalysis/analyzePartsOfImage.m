function analyzePartsOfImage(~,~,mainFig)


imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs'); 

% calculate spark frequency
imgDataXTfluoFN = imgData.imgDataXTfluoFN;
pxSzX = imgData.pxSzX;
pxSzT = imgData.pxSzT;
t = imgData.t;

if isfield(getappdata(mainFig),'sparkDetection') && isfield(getappdata(mainFig),'profileAnalysis')
    sparkDetection = getappdata(mainFig,'sparkDetection');
    detectedEvents = sparkDetection.detectedEvents;
    profileAnalysis = getappdata(mainFig,'profileAnalysis');
else
    errordlg('detect events and select profile!')
    return
end

% get size of part of image in s
w_t = str2double(hObjs.h_edit_imgSld_w.String);

% get start time of all events
s_E = zeros(size(detectedEvents,1),1);
for i=1:size(detectedEvents,1)   
    s_E(i,1) = detectedEvents(i).BoundingBox(1)*pxSzT;   
end

% get start time of selected repetitive events
s_rE = cell2mat(profileAnalysis.selectedROIs.eventsPeaks{1}(:,2));

t_idx = ceil(t./(w_t*1000));
N_parts = ceil(t(end)/(w_t*1000));

s = 1;
for i = 1:N_parts
    
    e = find(t_idx==i,1,'last');
    
    % calc area 
    imgArea = (size(imgDataXTfluoFN,1)*pxSzX)*(t(e)-t(s))/1000; % um*s
    % mask of events in this part of image
    mE = s_E>t(s) & s_E<=t(e);
    
    % calculate spark frequency
    sparkFreq(i,1) = sum(mE)*100/imgArea; % sparkFreq per 100um*s
    
    % calculate median of spark-to-spark delays
    m_RepSp = s_rE>t(s) & s_rE<=t(e);
    
    medParts(i,1) = median(diff(s_rE(m_RepSp)));
    
    s = e;
end

% plot results
hf = figure('Tag','parts of image','Units','normalized');
hf.Position = [0.25 0.5 0.5 0.5];
ax = axes('Parent',hf);
bar((1*w_t/2:w_t:numel(sparkFreq)*w_t),sparkFreq,0.9,'FaceColor','k','EdgeColor','k','FaceAlpha',0.25)
%line((1*w_t/2:w_t:numel(sparkFreq)*w_t),sparkFreq,'Linestyle','-','Marker','o','MarkerSize',10,'Color','k')
ax.YAxis.FontSize = 16;
ylabel('spark frequency (per 100 um*s)','FontSize',20)
ax.YAxis.Color = 'k';
ax.YLim = [0 1.05*max(sparkFreq)];
yyaxis right
bar((1*w_t/2:w_t:numel(medParts)*w_t),medParts,0.5,'FaceColor','r','EdgeColor','r','FaceAlpha',0.5)
%line((1*w_t/2:w_t:numel(medParts)*w_t),medParts,'Linestyle','-','Marker','o','MarkerSize',10,'Color','r')
ax.YAxis(2).FontSize = 16;
ax.XAxis.FontSize = 16;
ylabel('median of spark-to-spark delay (ms)','FontSize',20)
xlabel('t (s)','FontSize',20)
ax.YAxis(2).Color = 'r';
ax.YLim = [0 1.05*max(medParts)];

text(ax,min(ax.XLim),max(ax.YLim),{'width of part';sprintf('of image = %d s',w_t)},...
    'VerticalAlignment','top','FontSize',16)

% save data, spark frequencies and median of sp-to-sp delay
imgPartsAnalysis = struct('withOfPart_s',w_t,'sparkFreq_100um_s',sparkFreq,'spToSpDelay_med_ms',medParts);
profileAnalysis.selectedROIs.imgPartsAnalysis = imgPartsAnalysis;
setappdata(mainFig,'profileAnalysis',profileAnalysis)
% keyboard


end

