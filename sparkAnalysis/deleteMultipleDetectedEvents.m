function deleteMultipleDetectedEvents(hO,E)
% function to select and delete multiple detected events
hO.Parent.Pointer = 'cross';
switch E.Button
    case 1
        try
            % current point of pointer from axes
            cp = hO.CurrentPoint;
            % create rectangle
            xinit = cp(1,1);
            yinit = cp(1,2);
            r = rectangle(hO, ...
                'Position',[xinit yinit xinit-xinit yinit-yinit], ...
                'EdgeColor',[0 0 1], ...
                'LineWidth',3, ...
                'Interruptible','off');
            % get events rois data
            % get events data
            sparkDetection = getappdata(hO.Parent, 'sparkDetection');
            % save original color of events rois
            data.prevEvntsColor = cell2mat(arrayfun(@(x) get(x,'EdgeColor'), ...
                sparkDetection.detectedEventsRec, ...
                'UniformOutput',false));
            setappdata(r, 'data', data)
            % r.UserData.prevEvntsColor = cell2mat(arrayfun(@(x) get(x,'EdgeColor'), ...
            %     sparkDetection.detectedEventsRec, ...
            %     'UniformOutput',false));
            % setup mouse moving fcn and buttonup fcn
            hO.Parent.WindowButtonMotionFcn = {@wbmcb,hO,r,xinit,yinit,sparkDetection};
            hO.Parent.WindowButtonUpFcn = {@wbucb,hO,r};
        catch
            return
        end
    otherwise
        return
end


    % function executed when move with mouse pointer
    function wbmcb(hf, ~, ha, r, xinit, yinit, sparkDetection)
        % change color events to original one
        data = getappdata(r, 'data');
        if isfield(data, 'mask')
            for i = 1:numel(sparkDetection.detectedEventsRec)
                if ~data.mask(i)
                    set(sparkDetection.detectedEventsRec(i), ...
                        'EdgeColor',data.prevEvntsColor(i,:))
                    set(sparkDetection.detectedEventsMask(i), ...
                        'EdgeColor',data.prevEvntsColor(i,:), ...
                        'FaceColor',data.prevEvntsColor(i,:))
                end
            end
        end
        % get position of new point
        cp = ha.CurrentPoint;
        % update rectangle position
        r_w = abs(cp(1,1) - xinit);
        r_h = abs(cp(1,2) - yinit);
        if cp(1,1) < xinit
            r_x = cp(1,1);
        else
            r_x = xinit;
        end
        if cp(1,2) < yinit
            r_y = cp(1,2);
        else
            r_y = yinit;
        end
        r.Position = [r_x r_y r_w r_h];
        drawnow
        % find mask of events which fall in selection rectangle
        selectedEvntsMask = arrayfun(@(x) rectint(x.Position, r.Position), ...
            sparkDetection.detectedEventsRec) > 0;
        % change color of selected events to blue
        for i = 1:numel(sparkDetection.detectedEventsRec)
            if selectedEvntsMask(i)
                set(sparkDetection.detectedEventsRec(i), ...
                    'EdgeColor',[0 0 1])
                set(sparkDetection.detectedEventsMask(i), ...
                    'EdgeColor',[0 0 1], 'FaceColor',[0 0 1])
            end
        end
        drawnow
        % save mask of selected events to rectangle data
        data.mask = selectedEvntsMask;
        setappdata(r, 'data', data)
    end


    % stop function, when mouse button is released
    function wbucb(hf, ~, ha, r)
        
        % change back functions
        %hf.Pointer = 'arrow';
        hf.WindowButtonMotionFcn = '';
        hf.WindowButtonUpFcn = '';
        h_img = findobj(ha, 'Type', 'Image');
        h_img.PickableParts = 'none';
        % get selected events mask
        data = getappdata(r, 'data');
        selectedEvntsMask = data.mask;
        % delete selection rectangle
        delete(r)
        ha.ButtonDownFcn = '';
        sparkDetection = getappdata(hf, 'sparkDetection');
        % delete selected events
        % answer = questdlg( ...
        %     'Would you like to delete selected events?', ...
        % 	'Delete selected events.', ...
        % 	'YES','NO','NO');
        % switch answer
        %     case 'YES'
                % delete event from detection output
                newDetectedEvents = sparkDetection.detectedEvents;
                newDetectedEvents(selectedEvntsMask) = [];
                newTypeOfEvent = sparkDetection.typeOfEvent;
                newTypeOfEvent(selectedEvntsMask) = [];
                newMaskOfAcceptedSparks = sparkDetection.maskOfAcceptedSparks;
                newMaskOfAcceptedSparks(selectedEvntsMask) = [];
                sparkDetection.detectedEvents = newDetectedEvents;
                sparkDetection.typeOfEvent = newTypeOfEvent;
                sparkDetection.maskOfAcceptedSparks = newMaskOfAcceptedSparks;
                % save changes
                setappdata(hf, 'sparkDetection', sparkDetection)
                % update detected events
               eventsDetection(hO, [], hf, 'update')

        %     case 'NO'
        %         return
        % end
    end

end