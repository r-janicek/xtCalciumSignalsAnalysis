function showBiggerAxes(hO,E)

switch E.Button
    
    case 1 % create window
        
        % check if window exist
        hf = findall(0,'Tag','axesZoom');
        
        if isempty(hf)
   
            analysisFig = hO.Parent;
            hObjsA = getappdata(analysisFig,'hObjsA');
            
            hf = figure('Name','zoom of axes','Tag','axesZoom',...
                'units','normalized','outerposition',[0 0.5 1 0.5],...
                'CloseRequestFcn',{@closeAxesZoomFcn,analysisFig});
            
            % copy axes
            ax = copyobj(hO,hf);
            set(ax,'Position',[0.05 0.15 0.93 0.8])
            try
                mSz = get(0,'ScreenPixelsPerInch');
            catch 
                mSz = 50;
            end
            if hObjsA.check_expDPointSel.Value
                % find spline fit, get 
                hl_splFit = findall(ax,'Type','line','Tag','splineFit');
                [peakVal,peakPos] = max(hl_splFit.YData);

                sVal = hl_splFit.YData(1);
                sPos = 1;
                
                bs_ePos = find(hl_splFit.YData(1:peakPos) > ...
                    prctile(hl_splFit.YData, str2double(hObjsA.h_edit_bsSens.String) ), ...
                    1,'first');
                bs_eVal = hl_splFit.YData(bs_ePos);
                
                eVal = hl_splFit.YData(end);
                ePos = numel(hl_splFit.YData);

                % create points to select start/end for baseline fitting
                h_circ(1) = line(hl_splFit.XData(sPos),sVal,'Parent',ax,'Color','g',...
                    'LineStyle','none','Marker','.','MarkerSize',mSz,'PickableParts','all',...
                    'ButtonDownFcn',{@expDFitSelectionButtonDownFcn,analysisFig},...
                    'Tag','bsFitSelPoints');
                h_circ(2) = line(hl_splFit.XData(bs_ePos),bs_eVal,'Parent',ax,'Color','g',...
                    'LineStyle','none','Marker','.','MarkerSize',mSz,'PickableParts','all',...
                    'ButtonDownFcn',{@expDFitSelectionButtonDownFcn,analysisFig},...
                    'Tag','bsFitSelPoints');
                
                % create point to select peak
                h_circ(3) = line(hl_splFit.XData(peakPos),peakVal,'Parent',ax,'Color','m',...
                    'LineStyle','none','Marker','.','MarkerSize',mSz,'PickableParts','all',...
                    'ButtonDownFcn',{@expDFitSelectionButtonDownFcn,analysisFig},...
                    'Tag','peakFitSelPoints');
         
                % create points to select start/end for expD fitting
                h_circ(4) = line(hl_splFit.XData(peakPos),peakVal,'Parent',ax,'Color','k',...
                    'LineStyle','none','Marker','.','MarkerSize',mSz,'PickableParts','all',...
                    'ButtonDownFcn',{@expDFitSelectionButtonDownFcn,analysisFig},...
                    'Tag','expFitSelPoints');
                h_circ(5) = line(hl_splFit.XData(ePos),eVal,'Parent',ax,'Color','k',...
                    'LineStyle','none','Marker','.','MarkerSize',mSz,'PickableParts','all',...
                    'ButtonDownFcn',{@expDFitSelectionButtonDownFcn,analysisFig},...
                    'Tag','expFitSelPoints');
                
                % re-draw
                findAndAnalyzePeaks([],[],analysisFig)
                
                % move markers to top in axes
                uistack(h_circ,'top')
                
            end
            
        else
            % Figure exists so bring Figure to the focus
            figure(hf);
        end
        
        
    case 3 %delete window
        
        close(findobj(0,'Type','figure','Tag','axesZoom'))
        
end


end

