function saveImgForCOCOlikeDataset(~,~,mainFig)
% save processed image as OME-TIF
% save mask matrix in format uint8, where:
% 1 = spark, 2 = long-lasting spark, 3 = mini-wave, 4 = wave, 5 = transient

% change pointer to watch
set(mainFig,'Pointer','watch')
drawnow

imgData = getappdata(mainFig, 'imgData');
detectionResults = getappdata(mainFig, 'sparkDetection');

% pixelSize = ome.units.quantity.Length( java.lang.Double(imgData.pxSzX), ...
%    ome.units.UNITS.MICROMETER );
% pixelSizeT = ome.units.quantity.Time( java.lang.Double(imgData.pxSzT), ...
%    ome.units.UNITS.MILLISECOND );

% new name for part of image
[~,n,~] = fileparts(imgData.fileName);
% set up path to save processed image
% newImgPath = fullfile(imgData.filePath, [n, '_maskForML', '.ome.tiff']);
pathToSaveData = fullfile(imgData.filePath, [n, '_dataForML', '.mat']);

% create new metadata
% plane = uint16(imgData.imgDataXTfluoR);
% metadata = createMinimalOMEXMLMetadata(plane);
% metadata.setPixelsPhysicalSizeX(pixelSize, 0);
% metadata.setPixelsPhysicalSizeY(pixelSize, 0);
% metadata.setPlaneDeltaT(pixelSizeT,0,0);

% create mask of detected events
maskOfEvents = zeros(size(imgData.imgDataXTfluoR), 'uint8');
for i = 1:length( detectionResults.detectedEvents )
    if detectionResults.maskOfAcceptedSparks(i)
        maskOfEvents(detectionResults.detectedEvents(i).PixelIdxList) = ...
            detectionResults.typeOfEvent(i);
    end
end

% save image and mask as structure
savedImgData = struct('img', uint16(imgData.imgDataXTfluoR), ...
    'maskOfEvents', maskOfEvents, ...
    'pxSzX', imgData.pxSzX, ...
    'pxSzT', imgData.pxSzT);
save(pathToSaveData, 'savedImgData')

% put back pointer to arrow
set(mainFig,'Pointer','arrow')
drawnow  

end

