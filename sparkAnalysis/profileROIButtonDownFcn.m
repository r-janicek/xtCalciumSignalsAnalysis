function profileROIButtonDownFcn(h_O,event,mainFig)
% manage deleting detected events in cropped ROI or 
% delete and add detected points (local maxima) in profile
% load data
try
    hObjs = getappdata(mainFig,'hObjs');
    imgData = getappdata(mainFig,'imgData');
    profileAnalysis = getappdata(mainFig,'profileAnalysis');

    h_PeaksCirc = profileAnalysis.h_PeaksCirc;
    allCircPos = arrayfun(@(x) get(x,'XData'), h_PeaksCirc, ...
        'UniformOutput',1);
    posOfEvents = cell2mat(profileAnalysis.posOfEvents);
    smooth_span = str2double(get(hObjs.h_smooth_edit,'String'));
    bs_crit = str2double(get(hObjs.h_bsDet_edit,'String'));
    
    fitCoefSparkRise = profileAnalysis.fitCoefSparkRise;
    h_FitLine = profileAnalysis.h_FitLine;
       
    pxSzT = imgData.pxSzT;
    pxSzX = imgData.pxSzX;
    
    statEventsSpRec = profileAnalysis.statEventsSpRec;       
    croppedDataProfile = profileAnalysis.croppedDataWholeAreaRoi;
    
    ax_img_sparks = hObjs.ax_img_sparks;
    
    c = round(size(croppedDataProfile,1)/2);
    h_d = str2double(get(hObjs.h_edit_averageWidth,'String'));
    n_px = ceil(h_d/pxSzX);
    if mod(n_px,2)==0
        n_px = n_px - 1;
    end
    % create events mask for profile
    prof_t_evnts_m = false(size(croppedDataProfile,2),1);
    for i=1:length(statEventsSpRec)
        prof_t_evnts_m(statEventsSpRec(i).SubarrayIdx{2}, 1) = true;
    end
    % circles size
    cSz = profileAnalysis.peaksCircleSz;
    t = imgData.t;
    prof_t = get(findobj(hObjs.ax_prof, 'Type','Line', 'Color','k'),'YData');  
catch
    return
end


switch  h_O.Type
            
    case 'axes'
        
        switch event.Button
            
            case 1
                % create and mark new peak 
                
                point_pos = h_O.CurrentPoint;
       
                % find closet point in time profile
                [~,p] = min(abs(t - point_pos(1,1)));
                
                % look +- 10 ms and find peak if any
                r = ceil(20/pxSzT);
                
                % find peak
                try
                    [p_v,p_t] = findpeaks(prof_t(1,p-r:p+r),'SortStr','descend');
                    p_v = p_v(1);
                    p_t = p_t(1);
                    if isempty(p_v)
                        error('no peaks detected')
                    end
                catch
                    [p_v,p_t] = max(prof_t(1,p-r:p+r));
                end
                
                p_t = (p - r + p_t - 2)*pxSzT;
                
                % look if there is peak in this t position 
                try
                    [Lia,~] = ismember(p_t,allCircPos,'rows');
                catch
                    Lia = [];
                end
                
                if Lia
                    return
                end
                               
                h_circle = line(p_t, p_v, 'Parent',h_O, ...
                    'Color','r', 'LineStyle','none', ...
                    'Marker','.', 'MarkerSize',cSz, ...
                    'PickableParts','all',...
                    'ButtonDownFcn',{@profileROIButtonDownFcn,mainFig});               
                h_circ_n = [h_PeaksCirc; h_circle];
                
                posOfEvents = [posOfEvents;[p_v,p_t]];
                
                [posOfPeaks_n,index] = sortrows(posOfEvents,2);
                h_circ_n = h_circ_n(index);
                
                % do fit of events
                fitCoefSparkRise = [fitCoefSparkRise;[p_t-10 5 p_v 1]];
                fitCoefSparkRise = fitCoefSparkRise(index,:);
                pks = posOfEvents(:,1);
                pks = pks(index);
                locs = posOfEvents(:,2);
                locs = locs(index);
                delete(h_FitLine);
                               
                [h_line_n, detectedEventsMask_n, coef_n, ~,... 
                    startOfSpark,endOfSpark] = fitSparkRise(pxSzT, ...
                    t, prof_t, pks, locs, h_O, fitCoefSparkRise, ...
                    1e-3, 400, smooth_span, bs_crit, [], [], prof_t_evnts_m);
                              
                profileAnalysis.fitCoefSparkRise = coef_n;
                profileAnalysis.h_FitLine = h_line_n;
                profileAnalysis.detectedEventsMask = detectedEventsMask_n;
                profileAnalysis.startOfSpark = startOfSpark;
                profileAnalysis.endOfSpark = endOfSpark;
                profileAnalysis.h_PeaksCirc = h_circ_n;
                profileAnalysis.posOfEvents = num2cell(posOfPeaks_n);
                
                % save new data           
                setappdata(mainFig,'profileAnalysis',profileAnalysis);
                
            case 3
                return
        end
        
        
    case 'line'
        
        switch h_O.Marker
            
            case '.'               
                % delete selected peak 
                
                selectedCirclePos = h_O.XData;
                [Lia,Locb] = (ismember(selectedCirclePos,allCircPos,'rows'));
                
                switch event.Button
                    
                    case 1
                        return
                        
                    case 3
                        
                        if Lia
                            h_PeaksCirc(Locb) = [];
                            posOfEvents(Locb,:) = [];
                                                        
                            %do fit of events
                            fitCoefSparkRise(Locb,:) = [];
                            pks = posOfEvents(:,1);
                            locs = posOfEvents(:,2);
                          
                            delete(h_FitLine);
                            
                            % new fits of peaks                           
                            [h_line_n, detectedEventsMask_n, coef_n, ~, ...
                                startOfSpark, endOfSpark] = fitSparkRise( ...
                                pxSzT, t, prof_t, pks, locs, h_O.Parent, ...
                                fitCoefSparkRise, 1e-3, 400, ...
                                smooth_span, bs_crit, [], [], prof_t_evnts_m);
                            
                            profileAnalysis.fitCoefSparkRise = coef_n;
                            profileAnalysis.h_FitLine = h_line_n;
                            profileAnalysis.detectedEventsMask = detectedEventsMask_n;
                            profileAnalysis.startOfSpark = startOfSpark;
                            profileAnalysis.endOfSpark = endOfSpark;
                            profileAnalysis.h_PeaksCirc = h_PeaksCirc;
                            profileAnalysis.posOfEvents = num2cell(posOfEvents);
                            
                            % save new data
                            setappdata(mainFig,'profileAnalysis',profileAnalysis);

                            delete(h_O);
                            
                        end
                end
            

            case '+'
                
                switch event.Button
                    
                    case 1
                        return
                    
                    case 3
                        % delete selected detected event, in image
                        
                        % find index which object to delete
                        ind = find(h_O == [statEventsSpRec.centreLine]');                        
                        delete(statEventsSpRec(ind).eventRec)
                        statEventsSpRec(ind) = [];
                        delete(h_O)
                        profileAnalysis.statEventsSpRec = statEventsSpRec;
                        % save new data
                        setappdata(mainFig,'profileAnalysis',profileAnalysis);
                                             
                end
                
                
            case 'none'
                
                switch event.Button
                    
                    case 1

                        point_pos = h_O.Parent.CurrentPoint;
                        
                        % find closet point in time profile
                        [~,p] = min(abs(t - point_pos(1,1)));
                        
                        % look +- 10 ms and find peak if any
                        r = ceil(20/pxSzT);
                        
                        % find peak
                        try
                            [p_v,p_t] = findpeaks(prof_t(1,p-r:p+r),'SortStr','descend');
                            p_v = p_v(1);
                            p_t = p_t(1);
                            if isempty(p_v)
                                error('no peaks detected')
                            end
                        catch
                            [p_v,p_t] = max(prof_t(1,p-r:p+r));
                        end
                        
                        p_t = (p - r + p_t - 2)*pxSzT;
                        
                        % look if there is peak in this t position
                        try
                            [Lia,~] = ismember(p_t,allCircPos,'rows');
                        catch
                            Lia = [];
                        end
                        
                        if Lia
                            return
                        end
                        
                        h_circle = line(p_t,p_v,'Parent',h_O.Parent,'Color','r',...
                            'LineStyle','none','Marker','.','MarkerSize',cSz,'PickableParts','all',...
                            'ButtonDownFcn',{@profileROIButtonDownFcn,mainFig});
                        
                        h_circ_n = [h_PeaksCirc; h_circle];
                        
                        posOfEvents = [posOfEvents;[p_v,p_t]];
                        
                        [posOfPeaks_n,index] = sortrows(posOfEvents,2);
                        h_circ_n = h_circ_n(index);
                                                                     
                        % do fit of events
                        fitCoefSparkRise = [fitCoefSparkRise;[p_t-10 5 p_v 1]];
                        fitCoefSparkRise = fitCoefSparkRise(index,:);
                        pks = posOfEvents(:,1);
                        pks = pks(index);
                        locs = posOfEvents(:,2);
                        locs = locs(index);
                        delete(h_FitLine);                        
                        
                        [h_line_n,detectedEventsMask_n,coef_n,~,startOfSpark,endOfSpark] = ...
                            fitSparkRise(pxSzT,t,prof_t,pks,locs,h_O.Parent,fitCoefSparkRise,1e-3,400,smooth_span,bs_crit,[],[], prof_t_evnts_m);
                        
                        profileAnalysis.fitCoefSparkRise = coef_n;
                        profileAnalysis.h_FitLine = h_line_n;
                        profileAnalysis.detectedEventsMask = detectedEventsMask_n;
                        profileAnalysis.startOfSpark = startOfSpark;
                        profileAnalysis.endOfSpark = endOfSpark;
                        profileAnalysis.h_PeaksCirc = h_circ_n;
                        profileAnalysis.posOfEvents = num2cell(posOfPeaks_n);
                        
                        % save new data
                        setappdata(mainFig,'profileAnalysis',profileAnalysis);
                                               
                    case 3
                        return
                        
                end
                
        end
        
end

end


