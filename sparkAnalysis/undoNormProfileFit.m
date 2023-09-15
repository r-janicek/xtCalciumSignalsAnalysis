function undoNormProfileFit(~,~,fitFig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% get data
hObjsFit = getappdata(fitFig,'hObjsFit');
prof = getappdata(fitFig,'selectedProf');
profileAnalysis = getappdata(fitFig,'profileAnalysis');

% save normalized profile
if isfield(prof,'yN')   
    
    prof = rmfield(prof,'yN');      
    setappdata(fitFig,'selectedProf',prof);
    
    % setup fit window
    hObjsFit.popUpMenuBs.Enable = 'on';
    selectBsFitFcn(hObjsFit.popUpMenuBs,[],fitFig)
        
    hObjsFit.h_pb_fitBs.Enable = 'on';
    hObjsFit.h_pb_normProf.Enable = 'on';
    
    hObjsFit.maskButtonGroup.SelectedObject = hObjsFit.rbutton1;
    hObjsFit.rbutton1.Enable = 'on';
    hObjsFit.rbutton2.Enable = 'off';
    
    % plot profile and delete previous normalized profile and lines
    % delete previous fit
    delete(findobj(hObjsFit.ax_fit,'Tag','baselineFit'))
    delete(findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','Mask'))
    delete(findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','profile'))
    delete(findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','selectedPeaks'))
    delete(findobj(hObjsFit.ax_fit,'Type','Line'))
    
    % profile
    line(prof.t,prof.y,'Parent',hObjsFit.ax_fit,'Color','k','LineStyle','-','LineWidth',1,...
        'Tag','profile');
    set(hObjsFit.ax_fit,'XLim',[prof.t(1) prof.t(end)],'YLim',[min(prof.y)*0.95 max(prof.y)*1.05])
    
    % show selected peaks
    line([prof.eventsPeaks{:,2}],ones(size([prof.eventsPeaks{:,2}])).*hObjsFit.ax_fit.YLim(2),...
        'Parent',hObjsFit.ax_fit,'Color','r','LineStyle','none','LineWidth',1,...
        'Marker','.','MarkerSize',profileAnalysis.peaksCircleSz,'Tag','selectedPeaks');
    
    % adjust photolytic pulses if there are any
    PPulses = findobj(hObjsFit.ax_fit,'Tag','photolyticPulses');
    if ~isempty(PPulses)
        PPulses.YData = [ones(1,size(PPulses.YData,2)).*min(hObjsFit.ax_fit.YLim);
                        ones(1,size(PPulses.YData,2)).*min(hObjsFit.ax_fit.YLim);
                        ones(1,size(PPulses.YData,2)).*max(hObjsFit.ax_fit.YLim);
                        ones(1,size(PPulses.YData,2)).*max(hObjsFit.ax_fit.YLim)];
    
    end
     
    % do fitting of baseline
    fitBaseline([],[],fitFig)
       
end

end

