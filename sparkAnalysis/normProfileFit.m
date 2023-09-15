function normProfileFit(~,~,fitFig)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% get data
hObjsFit = getappdata(fitFig,'hObjsFit');
prof = getappdata(fitFig,'selectedProf');
imgData = getappdata(fitFig,'imgData');
profileAnalysis = getappdata(fitFig,'profileAnalysis');

F = prof.y(:);
F0 = prof.baselineFit(:);
b = imgData.blank;

% normalize selected profile (F-F0)/(F0-blank)
normProf = (F - F0) ./ (F0 - b);



% %%%%%%%%%%% 
% normProf = (F ./ F0) .* 2 .*112;

%.* mean( F0( find(prof.baselineM == 1,10,'last'))) ;
% 
% Fmin = mean(  normProf( round(119200/mean(diff(prof.t))  ) :round(120700/mean(diff(prof.t)) ) ) );
% 
% normProf = normProf .* ( 111/Fmin );
% 
% 
% % calcium concentration from calibration curve
% base = 111.87;
% maxF = 798.57;
% rate = 1.1514;
% xhalf = 0.32817;
% 
% 
% f_Ca = @ (base,maxF,rate,xhalf,F)  xhalf./( ( abs( (maxF-base)./(F-base) ) - 1 ).^(1/rate) );
% 
% % Ca in uM 
% normProf = f_Ca(base,maxF,rate,xhalf,normProf).*1000;
% 


% % 
% % normProf = normProf ./ mean(  normProf( round(119200/mean(diff(prof.t))  ) :round(120700/mean(diff(prof.t)) ) ) );
% % 
% % normProf = normProf - 1;
% 
% normProf = (F ./ F0) .* mean( F0( find(prof.baselineM == 1,10,'last')));
% 
% Fmin = mean(  normProf( round(119200/mean(diff(prof.t))  ) :round(120700/mean(diff(prof.t)) ) ) );
% Rm = 5;
  %keyboard
% 
% %%%%%%%%%%%



prof.yN = normProf;

% save normalized profile
setappdata(fitFig,'selectedProf',prof);

% setup fit window 
hObjsFit.popUpMenuBs.Enable = 'off';
hObjsFit.h_edit_paramFitBs1.Enable = 'off';
hObjsFit.h_edit_paramFitBs2.Enable = 'off';
hObjsFit.h_pb_fitBs.Enable = 'off';
hObjsFit.h_pb_normProf.Enable = 'off';

hObjsFit.maskButtonGroup.SelectedObject = hObjsFit.rbutton2;
hObjsFit.rbutton1.Enable = 'off';
hObjsFit.rbutton2.Enable = 'on';

% plot normalized profile and delete previous lines
% delete previous fit

delete(findobj(hObjsFit.ax_fit,'Tag','baselineFit'))
delete(findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','Mask'))
delete(findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','profile'))
delete(findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','selectedPeaks'))

% normalized profile 
line(prof.t,prof.yN,'Parent',hObjsFit.ax_fit,'Color','k','LineStyle','-','LineWidth',1,...
    'Tag','profile');
set(get(hObjsFit.ax_fit,'Ylabel'),'String',['profile (', (char(916)),'F/F0'])

% set ylimits
ymin = min(prof.yN);
if ymin<0
    ymin = ymin*1.05;
else
    ymin = ymin*0.95;
end

set(hObjsFit.ax_fit,'XLim',[prof.t(1) prof.t(end)],'YLim',[ymin max(prof.yN)*1.05])

% show selected peaks
line([prof.eventsPeaks{:,2}],ones(size([prof.eventsPeaks{:,2}])).*hObjsFit.ax_fit.YLim(2),...
    'Parent',hObjsFit.ax_fit,'Color','r','LineStyle','none','LineWidth',1,...
    'Marker','.','MarkerSize',profileAnalysis.peaksCircleSz,'Tag','selectedPeaks');

% show mask of events
hl_m = line(prof.t(prof.eventsM),prof.yN(prof.eventsM),'Parent',hObjsFit.ax_fit,'Color','g',...
                     'LineStyle','none','Marker','.','MarkerSize',20,'LineWidth',1,'Tag','eventsMask');
uistack(hl_m, 'bottom')

% adjust photolytic pulses if there are any
PPulses = findobj(hObjsFit.ax_fit,'Tag','photolyticPulses');
if ~isempty(PPulses)
    PPulses.YData = [ones(1,size(PPulses.YData,2)).*min(hObjsFit.ax_fit.YLim);
                     ones(1,size(PPulses.YData,2)).*min(hObjsFit.ax_fit.YLim);
                     ones(1,size(PPulses.YData,2)).*max(hObjsFit.ax_fit.YLim);
                     ones(1,size(PPulses.YData,2)).*max(hObjsFit.ax_fit.YLim)];
    
end

end

