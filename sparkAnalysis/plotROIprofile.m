function plotROIprofile(h_rect_prof_pos,mainFig)

set(mainFig,'Pointer','watch')
drawnow

% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
profileAnalysis = getappdata(mainFig,'profileAnalysis');
analysisType = getappdata(mainFig,'analysisType');

% get image data
switch analysisType
    
    case 'spark recovery ryanodine'
        img = imgData.imgDataXTfluoFN;
        
    case 'spark recovery photolysis'
        img = imgData.imgDataXTfluoF;
        
    otherwise
        img = imgData.imgDataXTfluoFN;
        
end

pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;
t = imgData.t;

ax_img_sparks = hObjs.ax_img_sparks;
%ax_img_sparks_2 = hObjs.ax_img_sparks_2;
crit = str2double(get(hObjs.h_edit_tresh_prof,'String'));
smooth_span = str2double(get(hObjs.h_smooth_edit,'String'));
bs_crit = str2double(get(hObjs.h_bsDet_edit,'String'));

ax_prof = hObjs.ax_prof;
h_rect_prof = profileAnalysis.h_rect_prof;

% get area of image for profile
% h_y_link = hObjs.hlink2;
% h_y_link.removetarget(ax_img_sparks)

h_rect_prof_pos = [floor(h_rect_prof_pos(1)) round(h_rect_prof_pos(2)) h_rect_prof_pos(3) h_rect_prof_pos(4)];

h_rect_prof_pos(1) = t(1);
h_rect_prof_pos(3) = t(end);

h_rect_prof.setPosition(h_rect_prof_pos)
%setappdata(main_fig,'h_rect_prof',h_rect_prof);

ROI_pos = [h_rect_prof_pos(1)/pxSzT h_rect_prof_pos(2) h_rect_prof_pos(3)/pxSzT+1 h_rect_prof_pos(4)];

ROI_pos(ROI_pos<0)=0;

[~,r_x(:,1)] = rect2ind(ROI_pos);

cropRoi = img(r_x,:); % crop from data for selected analysis 
cropRoiR = imgData.imgDataXTfluoR(r_x,:); % crop from raw data

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

image(cropRoi,'YData',[y_px(1) y_px(end)],'XData',[t(1) t(end)],...
    'CDataMapping','scaled','Parent',ax_img_sparks);
set(get(ax_img_sparks,'Ylabel'),'String','x (pixels)','FontWeight','bold', 'Tag','img_sparks')
set(ax_img_sparks,'XTick',[],'FontSize',14,'YDir','reverse','YAxisLocation','left')

% set and add interactive scale
sc_num = (diff(ax_img_sparks.YLim)*pxSzX)/2; % in um
hObjs.h_txt_scale.String = [sprintf('%0.2f',sc_num),' \mum'];
editScaleListener = addlistener(ax_img_sparks, 'YLim', 'PostSet', ...
    @(varargin)scaleLineChange(varargin,mainFig));
setappdata(mainFig,'editScaleListener',editScaleListener)

% show borders of area taken for average profile
line([t(1) t(end)],[c c],'Parent',ax_img_sparks,'LineWidth',3,'Color','k','LineStyle','-')
line([t(1) t(end)],[c-((n_px-1)/2) c-((n_px-1)/2)],'Parent',ax_img_sparks,'LineWidth',2,'Color','k','LineStyle',':')
line([t(1) t(end)],[c+((n_px-1)/2) c+((n_px-1)/2)],'Parent',ax_img_sparks,'LineWidth',2,'Color','k','LineStyle',':')


if isfield(imgData,'szImgs')
    if ~isempty(imgData.szImgs)
        % only for series of images, no analysis, no peaks and events
        % detection
        % take area for averaging, calculating time profile
        cropProf = cropRoi(c-((n_px-1)/2):c+((n_px-1)/2),:);
        prof_t(:,1) = mean(cropProf,1);
        
        % smooth profile with loess with defined duration in ms
        n_pts = round(smooth_span/pxSzT);
        prof_t_s = smooth(prof_t,n_pts/length(prof_t),'loess');
        
        % plot profile
        pp = plot(t,prof_t,t,prof_t_s,'Parent',ax_prof);
        set(pp(1),'Color','k')
        set(pp(2),'Color','b','LineStyle','-')
        if exist('oldXlim', 'var')
            xlim(ax_prof,[oldXlim(1) oldXlim(2)])
        else
            xlim(ax_prof,[min(t) max(t)])
        end
        ylim(ax_prof,[min(prof_t)*0.9 max(prof_t)*1.05])
        set(ax_prof,'FontSize',14)
        
        % set correct axis labels
        if isfield(imgData,'norm_flag') && imgData.norm_flag==1 && strcmp(analysisType,'spark recovery ryanodine')
            
            set(get(ax_prof,'Xlabel'),'String','t (px)','FontWeight','bold')
            set(get(ax_prof,'Ylabel'),'String',['fluorescence ','(',char(916),'F/F0)'],'FontWeight','bold')
            
        else
            set(get(ax_prof,'Xlabel'),'String','t (px)','FontWeight','bold')
            set(get(ax_prof,'Ylabel'),'String','fluorescence (F)','FontWeight','bold')
            
        end
        
        yh_patch = abs(diff(ax_prof.YLim))*0.10;
        
        szXimgCtrl = imgData.szImgs(imgData.m_ctrl,2);
        szXimgDrug = imgData.szImgs(~imgData.m_ctrl,2);
        
        cInc1 = 0.5/numel(szXimgCtrl);
        cInc2 = 0.6/numel(szXimgDrug);
        
        t_imgsCtrl = imgData.t_imgs(imgData.m_ctrl);
        t_imgsDrug = imgData.t_imgs(~imgData.m_ctrl);

        for i=1:numel(szXimgCtrl)
            
            ts = szXimgCtrl(i)*(i-1)+1;
            te = szXimgCtrl(i)*i;
            
            h_ctrlPatch = patch('XData',[ts te te ts],...
                'YData',[ax_prof.YLim(1) ax_prof.YLim(1) ax_prof.YLim(1)+yh_patch ax_prof.YLim(1)+yh_patch],...
                'EdgeColor','none','FaceColor','b','FaceAlpha',cInc1*i,'Parent',ax_prof);
            h_ctrlText = text(ts,ax_prof.YLim(1),sprintf('CTRL (img #%d, ts=%d s)',i,t_imgsCtrl(i)),'Color','w','FontSize',14,...
                'Parent',ax_prof,'VerticalAlignment','bottom','FontWeight','bold');
        end
        
        for i=1:numel(szXimgDrug)
            
            ts = szXimgDrug(i)*(i-1)+1+sum(szXimgCtrl);
            te = szXimgDrug(i)*i+sum(szXimgCtrl);
            
            h_drugPatch = patch('XData',[ts te te ts],...
                'YData',[ax_prof.YLim(1) ax_prof.YLim(1) ax_prof.YLim(1)+yh_patch ax_prof.YLim(1)+yh_patch],...
                'EdgeColor','none','FaceColor','r','FaceAlpha',cInc2*i,'Parent',ax_prof);
            h_drugText = text(ts,ax_prof.YLim(1),...
                sprintf('%s %d uM (img #%d, ts=%d s)',imgData.drugName,imgData.drugConc,i,t_imgsDrug(i)),'Color','k','FontSize',14,...
                'Parent',ax_prof,'VerticalAlignment','bottom','FontWeight','bold');
            
        end
    end
    
else
    % standard analysis
    % save cropped image data
    profileAnalysis.croppedDataWholeAreaRoi = cropRoi; % filtred and possibly normalized data
    profileAnalysis.croppedDataWholeAreaRoiR = cropRoiR; % raw data
    profileAnalysis.numOfPxAvrg = n_px;
    
    setappdata(mainFig,'profileAnalysis',profileAnalysis)
    
    %find individuals events
    try
        statEvents = findEvents(mainFig,cropRoiR);
    catch
        statEvents = [];
    end
    
    % remove events which bounding boxes are not crossing center lines
    try
        lines_px_r = (c-((n_px-1)/2):c+((n_px-1)/2));
        events_px_r = arrayfun(@(x) x.SubarrayIdx{1,1},statEvents,'UniformOutput',0);
        m_center = cellfun(@(x) any(ismember(lines_px_r,x)),events_px_r,'UniformOutput',1);
        statEvents(~m_center) = [];
    catch
        statEvents = [];
    end
    
    % find centre of the event
    pos_centre_x = arrayfun(@(x) x.WeightedCentroid(1,1),statEvents,'UniformOutput',1);
    pos_centre_y = round(arrayfun(@(x) x.WeightedCentroid(1,2),statEvents,'UniformOutput',1));
    
    %pos_max = num2cell([pos_centre_y,pos_centre_x],2);
    
    for i = 1:length(statEvents)
        
        event = cropRoi(statEvents(i).PixelIdxList);
        [~,p_m] = max(event);
        px = statEvents(i).PixelIdxList;
        
        col(i,1) = ceil(px(p_m)/size(cropRoi,1));
        row = rem(px(p_m),size(cropRoi,1));
        
    end
    
    try
        pos_max = num2cell([pos_centre_y,col],2);
    catch
        pos_max = [];
    end
    
    % show centre and boundary rectangle of detected events
    for i=1:length(statEvents)
        
        % bounding rectangle
        pos = statEvents(i).BoundingBox;
        eventRec = rectangle('Position',[pos(1)*pxSzT pos(2) pos(3)*pxSzT pos(4)],...
            'Parent',ax_img_sparks,'EdgeColor','r','LineStyle',':',...
            'LineWidth',2,'Tag',num2str(i));
        
        % centre line
        centre = pos_max{i};
        centreLine = line([centre(2)*pxSzT-pxSzT centre(2)*pxSzT-pxSzT],[centre(1) centre(1)],...
            'Parent',ax_img_sparks,'Color','r','LineStyle','none',...
            'Marker','+','MarkerSize',20,'LineWidth',3,'PickableParts','all',...
            'ButtonDownFcn',{@profileROIButtonDownFcn,mainFig});
        
        [statEvents(i).centreLine] = centreLine;
        [statEvents(i).eventRec] = eventRec;
        
    end
    
    % take area for averaging, calculating time profile
    cropProf = cropRoi(c-((n_px-1)/2):c+((n_px-1)/2),:);
    prof_t(:,1) = mean(cropProf,1);
    
    cropProfR = cropRoiR(c-((n_px-1)/2):c+((n_px-1)/2),:);
    % prof_t(:,1) = mean(cropProfR,1);
    
    profileAnalysis.croppedDataProfile = cropProf;
    profileAnalysis.croppedDataProfileR = cropProfR;
    
    % smooth profile with loess with defined duration in ms
    
    % %%%%%%%% cross validation, slow
    % nn = 500;
    % spans = linspace(0.0001,0.25,nn);
    % sse = zeros(size(spans));
    % cp = cvpartition(length(prof_t),'k',2);
    %
    % for i = 1:numel(spans)
    %
    %     f = @(train,test) norm(test(:,2) - mylowess(train,test(:,1),spans(i)))^2;
    %     sse(i) = sum( crossval(f,[t',prof_t],'partition',cp) );
    %     i
    % end
    %
    % [minsse,minj] = min(sse);
    % span = spans(minj);
    %
    % %%%%%%%%
    
    n_pts = round(smooth_span/pxSzT);
    prof_t_s = smooth(prof_t,n_pts/length(prof_t),'loess');
    
    % find peaks
    noise_std  = std(prof_t(prof_t < (mean(prof_t)+2*std(prof_t))));
    mph = mean(prof_t)+crit*noise_std; % minimum peak height
    mpw = 10; % minimum peak width
    mpd = mpw/2;% minimum peak distance
    mpp = 2*noise_std; % minimum peak prominence
    
    if max(prof_t) > mph
        
        [pks,locs] = findpeaks(prof_t,t,'MinPeakHeight',mph,...
            'MinPeakWidth',mpw,...
            'MinPeakDistance',mpd,...
            'MinPeakProminence',mpp);
    end
    
    % plot profile
    pp = plot(t,prof_t,t,prof_t_s,'Parent',ax_prof);
    set(pp(1),'Color','k')
    set(pp(2),'Color','b','LineStyle','-')
    
    if exist('oldXlim', 'var')
        xlim(ax_prof,[oldXlim(1) oldXlim(2)])
    else
        xlim(ax_prof,[min(t) max(t)])
    end
    ylim(ax_prof,[min(prof_t) max(prof_t)*1.1])
    set(ax_prof,'FontSize',14)
    
    set([ax_prof; pp], 'buttondownfcn', {@profileROIButtonDownFcn,mainFig})
    
    % scale circles
    scSz = get(0,'screensize');
    cSz = (min(scSz(3),scSz(4))/1440)*60;
    
    % plot detected peaks in profile
    if exist('pks','var')
        if ~isempty(pks)
            
            h_circ = zeros(length(pks),1);
            for i=1:length(pks)
                
                % plot detected peaks
                h_circ(i,1) = line(locs(i),pks(i),'Parent',ax_prof,'Color','r',...
                    'LineStyle','none','Marker','.','MarkerSize',cSz,'PickableParts','all',...
                    'ButtonDownFcn',{@profileROIButtonDownFcn,mainFig});
                
            end
         
            % fit sparks and plot results
            [h_line,detectedEventsMask,coef,~,startOfSpark,endOfSpark] = fitSparkRise(pxSzT,t,prof_t,pks,locs,ax_prof,[],1e-3,400,smooth_span,bs_crit,[],[]);
            
            
        end
        
    end
    
    
    % set correct axis labels
    if isfield(imgData,'norm_flag') && imgData.norm_flag==1 && strcmp(analysisType,'spark recovery ryanodine')
        
        set(get(ax_prof,'Xlabel'),'String','t (ms)','FontWeight','bold')
        set(get(ax_prof,'Ylabel'),'String',['fluorescence ' '(',char(916),'F/F0)'],'FontWeight','bold')
        
    else
        set(get(ax_prof,'Xlabel'),'String','t (ms)','FontWeight','bold')
        set(get(ax_prof,'Ylabel'),'String','fluorescence (F)','FontWeight','bold')
        
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
    
end


set(mainFig,'Pointer','arrow')
drawnow








%%%%%%%%%%%%%%%%%%%%
function ys=mylowess(xy,xs,span)
%MYLOWESS Lowess smoothing, preserving x values
%   YS=MYLOWESS(XY,XS) returns the smoothed version of the x/y data in the
%   two-column matrix XY, but evaluates the smooth at XS and returns the
%   smoothed values in YS.  Any values outside the range of XY are taken to
%   be equal to the closest values.

if nargin<3 || isempty(span)
    span = .3;
end

% Sort and get smoothed version of xy data
xy = sortrows(xy);
x1 = xy(:,1);
y1 = xy(:,2);
ys1 = smooth(x1,y1,span,'loess');

% Remove repeats so we can interpolate
mt = diff(x1)==0;
x1(mt)=[]; ys1(mt) = [];

% Interpolate to evaluate this at the xs values
ys = interp1(x1,ys1,xs,'linear',NaN);

% Some of the original points may have x values outside the range of the
% resampled data.  Those are now NaN because we could not interpolate them.
% Replace NaN by the closest smoothed value.  This amounts to extending the
% smooth curve using a horizontal line.
if any(isnan(ys))
    ys(xs<x1(1)) = ys1(1);
    ys(xs>x1(end)) = ys1(end);
end

end




end




