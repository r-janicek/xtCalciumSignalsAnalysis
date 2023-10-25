function mouseSetMaskFcn(hO,~)
try
    % setMask of events
    % change pointer
    hO.Parent.Pointer = 'cross';
    
    cp = hO.CurrentPoint;
    
    xinit = cp(1,1);
    yinit = cp(1,2);
    hl = line('XData',xinit,'YData',yinit,...
        'Marker','none' ,'color','g','LineWidth',5);
    
    hO.Parent.WindowButtonMotionFcn = {@wbmcb,hO,hl,xinit,yinit};
    hO.Parent.WindowButtonUpFcn = {@wbucb,hO,hl};
    
catch
    return
end


% function executed when move with mouse pointer
function wbmcb(hf,~,ha,hl,xinit,yinit)

cp = ha.CurrentPoint;
xdat = [xinit,cp(1,1)];
ydat = [yinit,cp(1,2)];
hl.XData = xdat;
hl.YData = ydat;
% set color of line, add = green, remove = red
if xinit < cp(1,1)
   hl.Color = 'g';   
else
   hl.Color = 'r';     
end

drawnow

end

% stop function, when mouse button is released
function wbucb(fig,~,ha,hl)

fig.Pointer = 'arrow';
fig.WindowButtonMotionFcn = '';
fig.WindowButtonUpFcn = '';

% set mask
hl_x = hl.XData;

switch getappdata(mainFig,'analysisType')
    case {'spark recovery ryanodine', 'spark detection'}
        hObjsFit = getappdata(fig,'hObjsFit');
        selectedProf = getappdata(fig,'selectedProf');
        switch hObjsFit.maskButtonGroup.SelectedObject.String

            case 'set mask of baseline'
                m = selectedProf.baselineM;
            case 'set mask of events'
                m = selectedProf.eventsM;
        end
        xt = selectedProf.t;
        if isfield(selectedProf,'yN')
            prof_t = selectedProf.yN;
        else
            prof_t = selectedProf.y;
        end

        if numel(hl_x)<2
            return
        end
        % position of start and end of line in points
        [~,s] = min(abs(xt - hl_x(1)));
        [~,e] = min(abs(xt - hl_x(2)));
        if s <= e
            % mark new area for event fitting
            m(s:e,1) = true;
        else
            % remove mark from selected area
            m(e:s,1) = false;
        end
        % delete selection line and marked events line
        delete(hl)
        eventsMaskLine = findobj(ha, 'Type','Line', ...
            '-regexp','Tag','Mask');
        delete(eventsMaskLine)
        % create new line
        try
            hl_m = line(xt(m),prof_t(m), 'Parent',ha, ...
                'Color','g', 'LineStyle','none', 'Marker','.', ...
                'MarkerSize',20, 'LineWidth',1, 'Tag','Mask');
            uistack(hl_m, 'bottom')
        catch
            keyboard
        end
        % save mask
        switch hObjsFit.maskButtonGroup.SelectedObject.String
            case 'set mask of baseline'
                selectedProf.baselineM = m;
                setappdata(fig,'selectedProf',selectedProf);
                % do fitting of baseline
                fitBaseline([],[],fig)
            case 'set mask of events'
                selectedProf.eventsM = m;
                setappdata(fig,'selectedProf',selectedProf);
        end

    case 'transients & waves'
        % get data
        imgData = getappdata(fig,'imgData');
        hObjs = getappdata(fig,'hObjs');
        m = imgData.baselineM;
        xt = imgData.t;
        profLineObj = findall(ha,'Tag','wholeCellProfile');
        prof_t = profLineObj.YData;
        if numel(hl_x)<2
            return
        end
        % position of start and end of line in points
        [~,s] = min(abs(xt - hl_x(1)));
        [~,e] = min(abs(xt - hl_x(2)));
        if s <= e
            % mark new area for event fitting
            m(s:e,1) = true;
        else
            % remove mark from selected area
            m(e:s,1) = false;
        end
        % delete selection line and marked events line
        delete(hl)
        eventsMaskLine = findobj(ha, 'Type','Line', ...
            '-regexp','Tag','Mask');
        delete(eventsMaskLine)
        % create new line
        try
            hl_m = line(xt(m),prof_t(m), 'Parent',ha, ...
                'Color','g', 'LineStyle','none', 'Marker','.', ...
                'MarkerSize',20, 'LineWidth',1, 'Tag','Mask');
            uistack(hl_m, 'bottom')
        catch
            keyboard
        end
        % do fitting of baseline
        type = hObjs.popUpMenuBaselineFcn.String{hObjs.popUpMenuBaselineFcn.Value};
        bsFit = fitBaseline(xt, prof_t, type, m, 1, ha);
        % save new baseline mask
        imgData.baselineM = m;
        imgData.bsFit = bsFit;
        setappdata(fig,'imgData',imgData);
end

end

end

