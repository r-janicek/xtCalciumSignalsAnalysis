function fitNum_slider(h_sld,~,fitFig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% get data
hObjsFit = getappdata(fitFig,'hObjsFit');
imgData = getappdata(fitFig,'imgData');
profileAnalysis = getappdata(fitFig,'profileAnalysis');
mainFig = getappdata(fitFig,'mainFig');

hObjsFit.h_txt_fitNum.String = sprintf('profile #%d',h_sld.Value);
hObjsFit.h_txt_fitNumSld.String = sprintf('profile #%d',h_sld.Value);

% check is selected profile was fitted, then only show it
if any(strcmp(profileAnalysis.selectedROIs.Properties.VariableNames,'wholeProfileFit'))
    
    if isstruct(profileAnalysis.selectedROIs.wholeProfileFit{h_sld.Value})
        isFitted = true;
    else
        isFitted = false;
    end
       
else
    isFitted = false;   
end


if isFitted
    
    % only show fit 
    fitData = profileAnalysis.selectedROIs.wholeProfileFit{h_sld.Value};
    
    % delete previous mask, profile, baseline fit and events fit
    delete(findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','Mask'))
    delete(findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','profile'))
    delete(findobj(hObjsFit.ax_fit,'Tag','baselineFit'))
    delete(findobj(hObjsFit.ax_fit,'Tag','eventsFit'))
    delete(findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','selectedPeaks'))
    
    % show normalized profile
    line(fitData.t,fitData.yN,'Parent',hObjsFit.ax_fit,'Color','k','LineStyle','-','LineWidth',1,...
        'Tag','profile');
    set(get(hObjsFit.ax_fit,'Ylabel'),'String',['profile (', (char(916)),'F/F0'])
    
    % set ylimits
    ymin = min(fitData.yN);
    if ymin<0
        ymin = ymin*1.05;
    else
        ymin = ymin*0.95;
    end
    
    set(hObjsFit.ax_fit,'XLim',[fitData.t(1) fitData.t(end)],'YLim',[ymin max(fitData.yN)*1.05])
    
    % show fit of profile
    line(fitData.t,fitData.profFit.wholeFit,'Parent',hObjsFit.ax_fit,'Color','r','LineStyle','-','LineWidth',2,...
        'Tag','eventsFit');
              
    % plot residuals
    delete(findobj(hObjsFit.ax_res,'Tag','residuals'))
    line(fitData.t,fitData.yN-fitData.profFit.wholeFit,'Parent',hObjsFit.ax_res,'Color','k','LineStyle','-','LineWidth',1,...
        'Tag','residuals');    
    % set ylimit
    ymin = min(fitData.yN-fitData.profFit.wholeFit);
    if ymin<0
        ymin = ymin*1.05;
    else
        ymin = ymin*0.95;
    end    
    set(hObjsFit.ax_res,'XLim',[fitData.t(1) fitData.t(end)],'YLim',[ymin max(fitData.yN-fitData.profFit.wholeFit)*1.05])
        
    % disable baseline fitting   
    hObjsFit.popUpMenuBs.Enable = 'off';
    hObjsFit.h_edit_paramFitBs1.Enable = 'off';
    hObjsFit.h_edit_paramFitBs2.Enable = 'off';
    hObjsFit.h_pb_fitBs.Enable = 'off';
    hObjsFit.h_pb_normProf.Enable = 'off';
    hObjsFit.h_pb_undoNormProf.Enable = 'off';
    
    hObjsFit.maskButtonGroup.SelectedObject = hObjsFit.rbutton2;
    hObjsFit.rbutton1.Enable = 'off';
    hObjsFit.rbutton2.Enable = 'on';
    
    % disable events fitting 
    hObjsFit.h_pb_fit.Enable = 'off';
    hObjsFit.h_pb_acceptFit.Enable = 'off';
        
else
    
    % get selected profile data
    prof = getSelectedProfileData(h_sld.Value,imgData,profileAnalysis,getappdata(mainFig,'analysisType'));
    
    % save data
    setappdata(fitFig,'selectedProf',prof) 
    
    % set up window
    % show selected profile
    delete(hObjsFit.ax_fit.Children)
    delete(hObjsFit.ax_res.Children)
    
    line(prof.t,prof.y,'Parent',hObjsFit.ax_fit,'Color','k','LineStyle','-','LineWidth',1,...
        'Tag','profile');
    set(hObjsFit.ax_fit,'XLim',[prof.t(1) prof.t(end)],'YLim',[min(prof.y)*0.95 max(prof.y)*1.05])
    set(get(hObjsFit.ax_fit,'Ylabel'),'String','profile (F)','FontWeight','bold')
    
    % show selected peaks
    line([prof.eventsPeaks{:,2}],ones(size([prof.eventsPeaks{:,2}])).*hObjsFit.ax_fit.YLim(2),...
        'Parent',hObjsFit.ax_fit,'Color','r','LineStyle','none','LineWidth',1,...
        'Marker','.','MarkerSize',profileAnalysis.peaksCircleSz,'Tag','selectedPeaks');

    
    % enable baseline fitting
    set(hObjsFit.popUpMenuBs,'Enable','on')   
    set(hObjsFit.h_edit_paramFitBs1,'Enable','on')
    set(hObjsFit.h_edit_paramFitBs2,'Enable','on')
    set(hObjsFit.h_pb_fitBs,'Enable','on')
    set(hObjsFit.h_pb_normProf,'Enable','on')
    hObjsFit.h_pb_undoNormProf.Enable = 'on';
    
    hObjsFit.maskButtonGroup.SelectedObject = hObjsFit.rbutton1;
    hObjsFit.rbutton1.Enable = 'on';
    hObjsFit.rbutton2.Enable = 'off';
    
    % enable events fitting 
    hObjsFit.h_pb_fit.Enable = 'on';
    hObjsFit.h_pb_acceptFit.Enable = 'on';
    
    % do fitting of baseline
    fitBaseline([],[],fitFig)
               
end

end
