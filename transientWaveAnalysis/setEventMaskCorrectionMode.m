function setEventMaskCorrectionMode(selectedItem, E, maskROI)

switch selectedItem.Text
    case 'add mask'
        maskROI.Color = [0 1 0];
        maskROI.Tag = 'addMask';

    case 'remove mask'
        maskROI.Color = [1 0 0];
        maskROI.Tag = 'removeMask';

end

end

