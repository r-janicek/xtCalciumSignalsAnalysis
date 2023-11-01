function fitNum_slider(h_sld, ~, analysisFig)

mainFig = findobj(0, 'Type','figure', 'Tag','mainFig');
switch getappdata(mainFig,'analysisType')
    case {'spark recovery ryanodine', 'spark detection'}
        % get data
        hObjsFit = getappdata(analysisFig,'hObjsFit');
        imgData = getappdata(analysisFig,'imgData');
        profileAnalysis = getappdata(analysisFig,'profileAnalysis');
        mainFig = getappdata(analysisFig,'mainFig');

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
            delete(findobj(hObjsFit.ax_fit, 'Type','Line', ...
                '-regexp','Tag','Mask'))
            delete(findobj(hObjsFit.ax_fit, 'Type','Line', ...
                '-regexp','Tag','profile'))
            delete(findobj(hObjsFit.ax_fit, 'Tag','baselineFit'))
            delete(findobj(hObjsFit.ax_fit, 'Tag','eventsFit'))
            delete(findobj(hObjsFit.ax_fit, 'Type','Line', ...
                '-regexp','Tag','selectedPeaks'))

            % show normalized profile
            line(fitData.t, fitData.yN, 'Parent',hObjsFit.ax_fit, ...
                'Color','k', 'LineStyle','-', 'LineWidth',1,...
                'Tag','profile');
            set(get(hObjsFit.ax_fit,'Ylabel'), ...
                'String',['profile (', (char(916)),'F/F0'])
            % set ylimits
            set(hObjsFit.ax_fit, ...
                'XLim',[fitData.t(1) fitData.t(end)], ...
                'YLim',getAxisLimits(fitData.yN, 5))
            % show fit of profile
            line(fitData.t, fitData.profFit.wholeFit, ...
                'Parent',hObjsFit.ax_fit, 'Color','r', ...
                'LineStyle','-', 'LineWidth',2,...
                'Tag','eventsFit');
            % plot residuals
            delete(findobj(hObjsFit.ax_res, 'Tag','residuals'))
            line(fitData.t, fitData.yN-fitData.profFit.wholeFit, ...
                'Parent',hObjsFit.ax_res, 'Color','k', ...
                'LineStyle','-', 'LineWidth',1,...
                'Tag','residuals');
            % set ylimit
            set(hObjsFit.ax_res, ...
                'XLim',[fitData.t(1) fitData.t(end)], ...
                'YLim',getAxisLimits(fitData.yN-fitData.profFit.wholeFit, 5))

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
            prof = getSelectedProfileData(h_sld.Value, ...
                imgData, profileAnalysis, ...
                getappdata(mainFig,'analysisType'));

            % save data
            setappdata(analysisFig, 'selectedProf', prof)

            % set up window
            % show selected profile
            delete(hObjsFit.ax_fit.Children)
            delete(hObjsFit.ax_res.Children)

            line(prof.t, prof.y, 'Parent',hObjsFit.ax_fit, ...
                'Color','k', 'LineStyle','-', 'LineWidth',1,...
                'Tag','profile');
            set(hObjsFit.ax_fit, ...
                'XLim',[prof.t(1) prof.t(end)], ...
                'YLim',getAxisLimits(prof.y, 5))
            set(get(hObjsFit.ax_fit,'Ylabel'), ...
                'String','profile (F)', 'FontWeight','bold')

            % show selected peaks
            line([prof.eventsPeaks{:,2}], ...
                ones(size([prof.eventsPeaks{:,2}])).*hObjsFit.ax_fit.YLim(2),...
                'Parent',hObjsFit.ax_fit, 'Color','r', ...
                'LineStyle','none', 'LineWidth',1, 'Marker','.', ...
                'MarkerSize',profileAnalysis.peaksCircleSz, ...
                'Tag','selectedPeaks');
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
            fitBaseline([],[],analysisFig)
        end


    case 'transients & waves'
        %% get data
        mainFig = getappdata(analysisFig,'mainFig');
        selectedROIs = getappdata(mainFig,'selectedROIs');
        hObjsA = getappdata(analysisFig,'hObjsA');
        imgData = getappdata(analysisFig,'imgData');

        % selected event number
        currEventNum = hObjsA.sld_fitNum.Value;

        % create selectedEvent structure
        selectedEvent.ROIdata = selectedROIs(currEventNum,:);
        selectedEvent.pxSzT = imgData.pxSzT;
        selectedEvent.pxSzX = imgData.pxSzX;
        roiName = selectedEvent.ROIdata.roiName{1};
        % save type of event
        if ~isempty(regexp(roiName,'wave','once'))
            selectedEvent.type = 'wave';
        elseif ~isempty(regexp(roiName,'caffeine','once'))
            selectedEvent.type = 'caffeine';
        else
            selectedEvent.type = 'transient';
        end
        % check if selected event was analyzed and accepted
        try
            isAnalyzed = selectedEvent.ROIdata.analysis{1}.accepted;
        catch
            isAnalyzed = false;
        end
        % save analysis
        if isAnalyzed
            selectedEvent.analysis = selectedEvent.ROIdata.analysis{1};
        end

        % save selected event
        setappdata(analysisFig, 'selectedEvent', selectedEvent)

        % set up window for fitting
        hObjsA.h_edit_paramFit1.String = ...
            num2str( ceil(max(selectedEvent.ROIdata.dataROIs.t)/50) );

        % set number of peaksfor a new event to 1
        hObjsA.h_edit_Npeaks.String = num2str(1);


        %% show filtered and normalized image of selected event
        % clear axes
        arrayfun(@(x) cla(x), findobj(analysisFig, 'Type','axes'))

        % show image
        image(selectedEvent.ROIdata.dataROIs.imgFN,...
            'YData',[min(selectedEvent.ROIdata.dataROIs.x) ...
                     max(selectedEvent.ROIdata.dataROIs.x)],...
            'XData',[min(selectedEvent.ROIdata.dataROIs.t) ...
                     max(selectedEvent.ROIdata.dataROIs.t)],...
            'CDataMapping','scaled', 'Parent',hObjsA.ax_orgImg);
        set(hObjsA.ax_orgImg, 'XTick',[], 'YGrid','off', 'FontSize',14)
        set(get(hObjsA.ax_orgImg,'Ylabel'), 'String','x (\mum)', ...
            'FontWeight','bold')
        set(get(hObjsA.ax_orgImg,'title'), ...
            'String','image filtered and normalized', 'FontWeight','bold')

        %% analyze image or show already analyzed image and accepted analysis
        if isAnalyzed
            % just show resutls, use analyze peaks function without fitting
            findAndAnalyzePeaks([], [], analysisFig)
        else
            % do whole analysis
            analyzeEvent([], [], analysisFig);
        end

        %% set up analysis window
        hObjsA.h_txt_fitNum.String = ...
            sprintf('analysis of event #%d:  %s', h_sld.Value, roiName);
        hObjsA.h_txt_fitNumSld.String = ...
            sprintf('event #%d;  %s', h_sld.Value, roiName);

        if isAnalyzed
            % set up window
            % disable fitting
            hObjsA.h_edit_paramFit1.Enable = 'off';
            hObjsA.h_edit_paramFit2.Enable = 'off';
            hObjsA.h_edit_Npeaks.Enable = 'off';
            hObjsA.check_signOfPeak.Enable = 'off';
            hObjsA.h_edit_bsSens.Enable = 'off';

            % disable mask correction
            hObjsA.h_edit_ellipseAxes1.Enable = 'off';
            hObjsA.h_edit_ellipseAxes2.Enable = 'off';
            hObjsA.h_edit_treshDet.Enable = 'off';
            hObjsA.h_edit_filtSzT.Enable = 'off';
            hObjsA.h_edit_filtSzX.Enable = 'off';
            hObjsA.h_edit_knotsSpan.Enable = 'off';
            hObjsA.h_pb_createROIcircle.Enable = 'off';

            % disable pushbuttons
            hObjsA.h_pb_acceptAnalysis.Enable = 'off';
            hObjsA.h_pb_refit.Enable = 'off';
            hObjsA.h_pb_detect.Enable = 'off';
            hObjsA.h_pb_deskew.Enable = 'off';
            hObjsA.h_pb_eventAnalysis.Enable = 'off';

        else
            % set up window
            % enable fitting
            hObjsA.h_edit_paramFit1.Enable = 'on';
            hObjsA.h_edit_paramFit2.Enable = 'on';
            hObjsA.h_edit_Npeaks.Enable = 'on';
            hObjsA.check_signOfPeak.Enable = 'on';
            hObjsA.h_edit_bsSens.Enable = 'on';

            % enable mask correction
            hObjsA.h_edit_ellipseAxes1.Enable = 'on';
            hObjsA.h_edit_ellipseAxes2.Enable = 'on';
            hObjsA.h_edit_treshDet.Enable = 'on';
            hObjsA.h_edit_filtSzT.Enable = 'on';
            hObjsA.h_edit_filtSzX.Enable = 'on';
            hObjsA.h_edit_knotsSpan.Enable = 'on';
            hObjsA.h_pb_createROIcircle.Enable = 'on';

            % enable pushbuttons
            hObjsA.h_pb_acceptAnalysis.Enable = 'on';
            hObjsA.h_pb_refit.Enable = 'on';
            hObjsA.h_pb_detect.Enable = 'on';
            hObjsA.h_pb_deskew.Enable = 'on';
            hObjsA.h_pb_eventAnalysis.Enable = 'on';

        end

end

end
