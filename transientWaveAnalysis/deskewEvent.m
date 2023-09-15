function deskewEvent(hO,E,analysisFig)

% data
hObjsA = getappdata(analysisFig,'hObjsA');
selectedEvent = getappdata(analysisFig,'selectedEvent');

%roiName = selectedROIs.roiName{indROI};
img = selectedEvent.ROIdata.dataROIs.imgFN;
eventImg = selectedEvent.detectedEvent.eventImg;
% pxSzT = selectedEvent.pxSzT;
% pxSzX = selectedEvent.pxSzX;
t = selectedEvent.ROIdata.dataROIs.t;
% x = selectedEvent.ROIdata.dataROIs.x;

% if there is active mask on detected event image use that
if isfield(selectedEvent.detectedEvent,'tempMask')
    
    % get image based on selected mask
    eventImg = zeros(size(img));
    eventImg(selectedEvent.detectedEvent.tempMask) = img(selectedEvent.detectedEvent.tempMask);
    
    % delete temporal axes
    delete(selectedEvent.detectedEvent.ax_maskROI)
   
    % save data 
    selectedEvent.detectedEvent.eventImg = eventImg;
    selectedEvent.detectedEvent.mask = selectedEvent.detectedEvent.tempMask;
    
    % show image of new detected event
    h_detImg = findall(hObjsA.ax_detEvent,'Type','image');
    h_detImg.CData = eventImg;

end
  
% get parameters
% % span for spline knots
k = str2double(hObjsA.h_edit_knotsSpan.String);    % in ms

% find minimum distance of detected wave from beggining of image,
% pixelwise, use 50% of peak of spline fit
% rearrange image to deskew wave
eventStartPx = zeros(size(eventImg,1),1);

nLines = size(eventImg,1);

% waitbar
hw = waitbar(0,'1','Name','deskewing...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(hw,'canceling',0);

for i=1:nLines
    % fit with spline
    yE = eventImg(i,:);
    y = img(i,:);
    % make sure y is positive, only for event start detection 
    y = y + abs(min(y(:)));
    ind_y_noZero = yE~=0;
    
    % setup knots every k ms
    nK = ceil( (max(t)-min(t))/k );
    splOrd = 3; % spline order, poly3
    % first guess of spline coefficients
    spl0  = spap2(nK, splOrd, t, y);
    % newknt for a possibly better knot distribution
    knots = newknt(spl0);
    spl0  = spap2(knots, splOrd, t, y);
    yf = fnval(spmak(spl0.knots,spl0.coefs),t);
    
    yf(~ind_y_noZero) = 0;
    bs = mean( y(1:find(ind_y_noZero~=0,1,'first')) );
    
    % check for gaps in detected signal

    % detect start of wave as 50 % of max amplitude
    [vMax,indMax] = max(yf);   
    ampl = vMax - bs;
    [~,eventStartPx(i)] = min( abs( yf(1:indMax) - (0.5*ampl+bs) ) );

%     if eventStartPx(i)==1
%         keyboard
%     end

    % Update waitbar and message
    waitbar(i/(2*nLines),hw,sprintf('%d %%',round((i/(2*nLines))*100)))
    
%     figure
%     plot(t,y)
%     hold on
%     plot(t(ind_y_noZero), y(ind_y_noZero),'og')
%     plot(t,yf,'r')
%     plot(t(eventStartPx(i)),yf(eventStartPx(i)),'*b')
     
end

minDist = min(eventStartPx);
deskewedEvent = nan(size(eventImg));
for i=1:size(eventImg,1)
    d = eventStartPx(i) - minDist + 1;
    deskewedEvent(i,:) = [img(i,d+1:end),img(i,1:d)];
    
    % Update waitbar and message
    waitbar((nLines+i)/(2*nLines),hw,sprintf('%d %%',round(((nLines+i)/(2*nLines))*100)))
end


% delete waitbar
delete(hw)


% save
selectedEvent.deskewdEvent.eventStartPx = eventStartPx;
selectedEvent.deskewdEvent.deskewedImg = deskewedEvent;
selectedEvent.deskewdEvent.deskewedImg_tProf = mean(deskewedEvent,1);

setappdata(analysisFig,'selectedEvent',selectedEvent)

end

