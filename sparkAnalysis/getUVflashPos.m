function flashPos = getUVflashPos(hO,E,mainFig)
% find position of UV flash

imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');
if isempty(imgData)
    return
end

if isempty(hO)
    hO.Style = 'open image fnc';
end

% delete previous lines
delete(findall(hObjs.ax_img, 'Tag','UVflashPos'))

if hObjs.h_check_UVflash.Value
    switch hO.Style
        case {'checkbox', 'open image fnc'}
            % get UV flash position in ms
            if ~isempty(imgData.imgDataXTtrans)
                [~,flashPos] = max( ...
                    gradient( ...
                    smooth(mean(imgData.imgDataXTtrans,1), 5) ) );
            else
                [~,flashPos] = max( ...
                    gradient( ...
                    smooth(mean(imgData.imgDataXTfluoFN,1), 5) ) );
            end
            flashPos = round(flashPos*imgData.pxSzT); % in ms
            hObjs.h_edit_UVflash.String = num2str(flashPos);

        case 'edit'
            flashPos = str2double(hObjs.h_edit_UVflash.String);
    end

    % show UV flash position as vertical magenta line in line scan
    line(hObjs.ax_img, [flashPos,flashPos], hObjs.ax_img.YLim, ...
        'Color','m', 'LineWidth',2, 'Tag','UVflashPos')

    % add UV flash position also to notes
    % check if there is any
    ind = find( ...
        contains(hObjs.h_table_notes.Data, 'UVflashPos'), 1, 'first');
    if ~isempty(ind)
        % rewrite
        hObjs.h_table_notes.Data(ind) = ...
            {sprintf('UVflashPos = %0.2f (ms)', flashPos)};
    else
        % add at the end of notes, leave one empty before
        ind_last_nonEmpty = find( ...
            cellfun(@(x) ~isempty(x), hObjs.h_table_notes.Data), 1, 'last');
        hObjs.h_table_notes.Data = [...
            hObjs.h_table_notes.Data(1:ind_last_nonEmpty); {''}; ...
            {sprintf('UVflashPos = %0.2f (ms)', flashPos)}];
    end

else
    % remove UVflasPos from notes and set edit field to nan
    ind = find(contains(hObjs.h_table_notes.Data, 'UVflashPos'),1);
    if ~isempty(ind)
        hObjs.h_table_notes.Data(ind) = {''};
    end
    hObjs.h_edit_UVflash.String = 'nan';
end







