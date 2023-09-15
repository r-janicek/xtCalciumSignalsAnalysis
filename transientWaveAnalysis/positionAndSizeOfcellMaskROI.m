function positionAndSizeOfcellMaskROI(hFig, event, maskROI)

% get position of cursor in figure relatively to axes with mask
selectedEvent = getappdata(hFig,'selectedEvent');
ax = selectedEvent.detectedEvent.ax_maskROI;
if ~isgraphics(ax)
    return
end
curPos = get(ax, 'CurrentPoint');

% check if cursor is above the axes
xCurPos = curPos(1,1);
yCurPos = curPos(1,2);
if (xCurPos<=ax.XLim(2)) && (xCurPos>=ax.XLim(1)) && ...
        (yCurPos<=ax.YLim(2)) && (yCurPos>=ax.YLim(1))
    % update roi position or size
    switch event.EventName
        case 'WindowMouseMotion'
            maskROI.Position = [...
                xCurPos-maskROI.Position(3)/2,...
                yCurPos-maskROI.Position(4)/2,...
                maskROI.Position(3),...
                maskROI.Position(4)
                ];

        case 'WindowScrollWheel'
            if ax.XLim(2)>ax.YLim(2)
                x_mag = 1;
                t_mag = ceil(ax.XLim(2)/ax.YLim(2));
            else
                x_mag = ceil(ax.YLim(2)/ax.XLim(2));
                t_mag = 1;
            end
            szROI_x = maskROI.Position(4)+x_mag*event.VerticalScrollCount;
            szROI_t = maskROI.Position(3)+t_mag*event.VerticalScrollCount;
            if szROI_x<3
                szROI_x=3;
            elseif szROI_x>ax.YLim(2)
                szROI_x = ax.YLim(2);
            end
            if szROI_t<3
                szROI_t=3;
            elseif szROI_x>ax.XLim(2)
                szROI_t = ax.XLim(2);
            end
            % set ROI size
            maskROI.Position = [...
                xCurPos-maskROI.Position(3)/2,...
                yCurPos-maskROI.Position(4)/2,...
                szROI_t,...
                szROI_x
                ];

    end
    % add pause, not to overload
    pause(0.0001)
end





end