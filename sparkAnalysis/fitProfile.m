function fitProfile(~,~,fitFig)

set(fitFig,'Pointer','watch')
drawnow

% get data
hObjsFit = getappdata(fitFig,'hObjsFit');
prof = getappdata(fitFig,'selectedProf');
imgData = getappdata(fitFig,'imgData');

% get parameters for events fitting
pos_peak = round([prof.eventsPeaks{:,2}]./imgData.pxSzT);
t = prof.t;
y = prof.yN;
baselineModel = 'poly1';
fitFun = hObjsFit.popUpMenuEventsFcn.String{hObjsFit.popUpMenuEventsFcn.Value};
t_ups = (t(1):imgData.pxSzT/10:t(end));
recalculate = 0;
useWeights = hObjsFit.check_fitWeights.Value;
mEvents = prof.eventsM;
mBaseline = prof.baselineM;
setappdata(fitFig, 'stopFiting', false)

% spark rise coeff = [t0,tR,A,y0]
prevFitCoefSpRise = prof.prevFitCoefSpRise;

% delete previous fit
delete(findobj(hObjsFit.ax_fit,'Tag','eventsFit'))
delete(findobj(hObjsFit.ax_fit,'Tag','finalPartialFit'))

% do fit 
profFit = multipleEventsFitFcnNormBs(pos_peak, t, y,baselineModel,...
     fitFun, t_ups, recalculate, mEvents, mBaseline, hObjsFit.ax_fit,...
     useWeights, prevFitCoefSpRise);

% profFit = multipleEventsFitFcnNormBs(pos_peak,t,y,baselineModel,...
%     fitFun,t_ups,recalculate,mEvents,mBaseline,hObjsFit.ax_fit,useWeights);

% delete partial fits
delete(findobj(hObjsFit.ax_fit,'Tag','finalPartialFit'))

% plot results 
% new whole fit
line(t_ups,profFit.t_ups.wholeFit, ...
    'Parent',hObjsFit.ax_fit, 'Color','r', 'LineStyle','-', ...
    'LineWidth',2, 'Tag','eventsFit');
% partial fits
% too slow just show only part of data, do not show 0(baseline), 
% visible off is working
allIndEvents = cell2mat(profFit.t_ups.individualEventsFits(2:end,:));
allIndEvents(allIndEvents<0.001) = nan;
line(t_ups, allIndEvents, 'Parent',hObjsFit.ax_fit, ...
    'Color','b', 'LineStyle','-', 'LineWidth',1,...
    'Tag','eventsFit','Visible','on');         
    
% delete previous mask
hl_m_Old = findobj(hObjsFit.ax_fit, ...
    'Type','Line', '-regexp','Tag','Mask');
delete(hl_m_Old)
% % new mask
hl_m = line(t(mEvents), y(mEvents), 'Parent',hObjsFit.ax_fit, ...
    'Color','g', 'LineStyle','none', 'Marker','.', ...
    'MarkerSize',20, 'LineWidth',1, 'Tag','eventsMask');
uistack(hl_m, 'bottom')

% plot residuals
delete(findobj(hObjsFit.ax_res, 'Tag','residuals'))
line(t, y-profFit.wholeFit, 'Parent',hObjsFit.ax_res, ...
    'Color','k', 'LineStyle','-', 'LineWidth',1, 'Tag','residuals');
% set ylimit
set(hObjsFit.ax_res, 'XLim',[t(1) t(end)], ...
    'YLim',getAxisLimits(y-profFit.wholeFit, 5))

% save data
prof.profFit = profFit;
setappdata(fitFig, 'selectedProf', prof);

set(fitFig, 'Pointer', 'arrow')
drawnow

end

