function analyzeEvent(hO,~,analysisFig)

set(analysisFig,'Pointer','watch')
drawnow

selectedEvent = getappdata(analysisFig,'selectedEvent');
% remove previously set figure functions
analysisFig.WindowScrollWheelFcn = [];
analysisFig.WindowButtonMotionFcn = [];

% data
hObjsA = getappdata(analysisFig,'hObjsA');
roiName = selectedEvent.ROIdata.roiName{1};

img = selectedEvent.ROIdata.dataROIs.imgFN;
t = selectedEvent.ROIdata.dataROIs.t;
x = selectedEvent.ROIdata.dataROIs.x;
pxSzT = selectedEvent.pxSzT;
pxSzX = selectedEvent.pxSzX;

switch selectedEvent.type
    
    case 'wave'
    %% %%%%%%%%%%%%%%%%%%%%%% WAVE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % detect event
    detectEvent(hO,[],analysisFig);
    
    % if it is only changing mask of detected image, basicaly calling function
    % from correct mask event - parameters panel
    if ~isempty(hO)
        if strcmp(hO.Parent.Title,'parameters')
            set(analysisFig,'Pointer','arrow')
            drawnow
            return
        end
    end
    
    % deskew event
    deskewEvent(hO,[],analysisFig);
    
    selectedEvent = getappdata(analysisFig,'selectedEvent');
    waveImgDeskewed = selectedEvent.deskewdEvent.deskewedImg;
    waveStartPx = selectedEvent.deskewdEvent.eventStartPx;
    waveImg = selectedEvent.detectedEvent.eventImg;
    
    % t-profile
    deskewedWaveProf = selectedEvent.deskewdEvent.deskewedImg_tProf;
    
    
    % get exp decay tau based on different amplitude regions in wave
    [~,maxTprof]= max(deskewedWaveProf);

    % take ~10 ms wide stripe, odd number of pixels
    nPxT = round(10/pxSzT);
    if mod(nPxT,2)~=0
        nPxT = nPxT+1;
    end
    if (maxTprof-nPxT/2)<1
        waveProfX = mean(waveImgDeskewed(:,1:maxTprof+nPxT/2),2);
    elseif (maxTprof+nPxT/2)>size(waveImgDeskewed,2)
        waveProfX = mean(waveImgDeskewed(:,maxTprof-nPxT/2:end),2);
    else
        waveProfX = mean(waveImgDeskewed(:,maxTprof-nPxT/2:maxTprof+nPxT/2),2);
    end
    % fit with spline
    nK = ceil( (max(x)-min(x))/1 ); % fit every 1 um
    splOrd = 3; % spline order, poly3
    spl0  = spap2(nK, splOrd, x, waveProfX);  % first guess of spline coefficients
    % newknt for a possibly better knot distribution
    knots = newknt(spl0);
    
    spl0  = spap2(knots, splOrd, x, waveProfX);
    yf_x = fnval(spmak(spl0.knots,spl0.coefs),x);
    
    % find separation between regions with different amplitudes
    [~,minLocs] = findpeaks(-yf_x);
    
    % calculate exp decay time constant for different regions 
    % exp decay function 
    funExpFit = @(p,x,y) (p(1)*exp(-(x-p(2))/p(3))+p(4)) - y;
    funExp = @(p,x) p(1)*exp(-(x-p(2))/p(3))+p(4);
    optFit = optimoptions('lsqnonlin','Display','off');
    
    tauD_regionWise = num2cell(zeros(numel(minLocs)+1,3));
    nRegions = numel(minLocs)+1;
    for i=1:nRegions
        if i==1
            tProfReg = mean(waveImgDeskewed(1:minLocs(i),:),1);
            regPx = [1 minLocs(i)];% region pixels
            regU =  regPx.*pxSzX;% region units, um
        elseif i==nRegions
            tProfReg = mean(waveImgDeskewed(minLocs(i-1):end,:),1);
            regPx = [minLocs(i-1) size(waveImgDeskewed,1)];% region pixels
            regU =  regPx.*pxSzX;% region units, um
        else
            tProfReg = mean(waveImgDeskewed(minLocs(i-1):minLocs(i),:),1);
            regPx = [minLocs(i-1) minLocs(i)];% region pixels
            regU =  regPx.*pxSzX;% region units, um
        end
        
        % fit with spline
        nK = ceil( (max(t)-min(t))/20 ); % fit every 20 ms
        splOrd = 3; % spline order, poly3
        spl0  = spap2(nK, splOrd, t, tProfReg);  % first guess of spline coefficients
        % newknt for a possibly better knot distribution
        knots = newknt(spl0);
        spl0  = spap2(knots, splOrd, t, tProfReg);
        yf = fnval(spmak(spl0.knots,spl0.coefs),t);
        
        % find tr*max in decay part
        [vMax,indMax] = max(yf);
        
        tr = 0.80;
        [~,expDstartFit] = min( abs( yf(indMax:end) - tr*vMax ) );
    
        % exp decay fit 
        t_exp = t(indMax+expDstartFit-1:end);
        y_exp = tProfReg(indMax+expDstartFit-1:end);
        try
            coefExp = lsqnonlin(@(p)funExpFit(p,t_exp,y_exp),...
                [vMax,t_exp(1),20,0],[0 0 0 -inf],[inf inf inf inf],optFit);
        catch
            coefExp = nan(1,4);
        end
        % tauD in ms; position of region in pixels and in um
        tauD_regionWise(i,1) = {coefExp(3)};
        tauD_regionWise(i,2) = {regPx};
        tauD_regionWise(i,3) = {regU};
         
%         figure
%         plot(t,tProfReg,'ok')
%         hold on
%         plot(t,yf,'r')    
%         plot(t_exp,y_exp,'b')
%         plot(t_exp,funExp(coefExp,t_exp),'m','LineWidth',2)
%        
    end
    
    
    % get wave speed
    waveStartT = t(waveStartPx)./1000; % in s 
    waveStartT = waveStartT(:);
    x = x(:);
    
    % fit with spline with order 2 (piece-wise line)
    % first guess of spline coefficients
    spl = spap2(ceil(max(x)/5), 2, x, waveStartT);
    % newknt for a possibly better knot distribution
    knots = newknt(spl);
    spl  = spap2(knots, 2, x, waveStartT);
    splFit = fnval(spl,x);
    
    % 1. derivations of spline
    dx = gradient(x);
    dt = gradient(splFit);
    
    % add time dt into calculation
%      splFit_der = fnder(spl,1);  
%      y_prime = fnval(splFit_der,x); 
%     
    waveSpeedPx = abs(dx./dt); % um/s
    
    
%     keyboard
%     
%     figure
%     plot(waveSpeedPx,x,'ok')
%     hold on
%     plot(splFit,x,'r')
    
    % check for outliers
    mOutlWsp = isoutlier(waveSpeedPx,'median');
    % calculate wave speed without outliers
    waveSpeed = median(waveSpeedPx(~mOutlWsp));
    
    % show line to calculate wave speed
    delete(findobj(hObjsA.ax_orgImg,'Tag','fitDeskew'))

    br = find( [0;diff(mOutlWsp)] ~= 0 );
    for i = 1:numel(br)+1
        if i==1
            if i==numel(br)+1
                % only one segment, br is empty
                xp = x;
                yp = splFit.*1000; % in ms
                mOutlP = mOutlWsp;
            else
                % segments with same characteristics
                xp = x(1:br(i)-1);
                yp = splFit(1:br(i)-1).*1000; % in ms
                mOutlP = mOutlWsp(1:br(i)-1);
            end
            
        elseif i==numel(br)+1
            xp = x(br(i-1)-1:end);
            yp = splFit(br(i-1)-1:end).*1000; % in ms
            mOutlP = mOutlWsp(br(i-1)-1:end);
            
        else
            xp = x(br(i-1):br(i)-1);
            yp = splFit(br(i-1):br(i)-1).*1000; % in ms
            mOutlP = mOutlWsp(br(i-1):br(i)-1);
        end
        
        if mode(mOutlP) == false
            lColor = [0 1 0 0.75];
        else
            lColor = [1 0 0 0.75];
        end
        
        line(yp,xp,...
                'Parent',hObjsA.ax_orgImg,'Color',lColor,...
                'LineStyle','-','LineWidth',3,...
                'Tag','fitDeskew');
                      
    end
    
%     line(splFit(~mOutlWsp)*1000,x(~mOutlWsp),...
%         'Parent',hObjsA.ax_orgImg,'Color',[0 1 0 0.75],'LineStyle','-','LineWidth',3,...
%         'Tag','fitDeskew');
%     
%     line(splFit(mOutlWsp)*1000,x(mOutlWsp),...
%         'Parent',hObjsA.ax_orgImg,'Color',[1 0 0 0.75],'LineStyle','-','LineWidth',3,...
%         'Tag','fitDeskew');
     
    fitWaveFront.waveStartT = waveStartT; % s
    fitWaveFront.x = x; % um
    fitWaveFront.splFit = splFit;
    fitWaveFront.mOutlWsp = mOutlWsp;
    
    
    % clear axes
    delete(hObjsA.ax_detEvent.Children)
    delete(hObjsA.ax_deskewed.Children)
    
    % show analysis
    image(waveImg,'YData',[min(x) max(x)],...
        'XData',[min(t) max(t)],...
        'CDataMapping','scaled','Parent',hObjsA.ax_detEvent);
    set(hObjsA.ax_detEvent,'YGrid','off','FontSize',14)
    set(get(hObjsA.ax_detEvent,'Ylabel'),'String','x (\mum)','FontWeight','bold')
    set(get(hObjsA.ax_detEvent,'Xlabel'),'String','t (ms)','FontWeight','bold')
    set(get(hObjsA.ax_detEvent,'title'),'String','detected event','FontWeight','bold')
    
    image(waveImgDeskewed,'YData',[min(x) max(x)],...
        'XData',[min(t) max(t)],...
        'CDataMapping','scaled','Parent',hObjsA.ax_deskewed);
    set(hObjsA.ax_deskewed,'XTick',[],'YGrid','off','FontSize',14)
    set(get(hObjsA.ax_deskewed,'Ylabel'),'String','x (\mum)','FontWeight','bold')
    set(get(hObjsA.ax_deskewed,'title'),'String','deskewed event','FontWeight','bold')
    
    % set x axis
    hObjsA.ax_detEvent.XLim = [min(t) max(t)];

%     % show different amplitudes regions 
%     tauD_regionWise
%     
%     regLinesY = unique(reshape(regLinesY',[],1));
%   
%     line(hObjsA.ax_deskewed.XLim, [regLinesY regLinesY],...
%         'Parent',hObjsA.ax_deskewed,'Color','k','LineStyle',...
%         '-','LineWidth',1,'Tag','amplRegions');
%     
%     line(yf_x,x,...
%         'Parent',hObjsA.ax_deskewed,'Color','r','LineStyle',...
%         '-','LineWidth',1,'Tag','amplRegions');
%   
%     delete( findall(0,'Type','Line','Tag','amplRegions')  )
    
    % show deskewed profile
    line(t,deskewedWaveProf,...
        'Parent',hObjsA.ax_deskewedProf,'Color',[0 0 0 0.75],'LineStyle','-','LineWidth',1,...
        'Tag','profile');
    
    % save data
    waveAnalysis = struct('waveImg',waveImg,'waveImgDeskewed',waveImgDeskewed,...
                          'tauD_regionWise',{tauD_regionWise},'waveSpeed',waveSpeed,...
                          'fitOfWaveFront',{fitWaveFront},...
                          'deskewedWaveProf',deskewedWaveProf);
    
    selectedEvent.analysis.waveAnalysis = waveAnalysis;
    selectedEvent.analysis.transientAnalysis = [];
    selectedEvent.analysis.caffeineAnalysis = [];
    selectedEvent.analysis.accepted = false;
    
    setappdata(analysisFig,'selectedEvent',selectedEvent)
    
    % fit t-profile of deskewed wave with spline and with exponential rise
    % to get t0
    findAndAnalyzePeaks([],[],analysisFig)
    
    % set up zoom
    hObjsA.ax_deskewedProf.ButtonDownFcn = @showBiggerAxes;
    
    
      
case 'caffeine' 
    %% %%%%%%%%%%%%%%%%%%%% CAFFEINE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % detect event
    detectEvent(hO,[],analysisFig);
    
    % if it is only changing mask of detected image, basicaly calling function
    % from correct mask event - parameters panel
    if ~isempty(hO)
        if strcmp(hO.Parent.Title,'parameters')
            set(analysisFig,'Pointer','arrow')
            drawnow
            return
        end
    end
    
    % deskew event
    deskewEvent(hO,[],analysisFig);
    
    selectedEvent = getappdata(analysisFig,'selectedEvent');
    waveImgDeskewed = selectedEvent.deskewdEvent.deskewedImg;
    waveStartPx = selectedEvent.deskewdEvent.eventStartPx;
    waveImg = selectedEvent.detectedEvent.eventImg;
    
    % t-profile
    deskewedWaveProf = selectedEvent.deskewdEvent.deskewedImg_tProf;
    
    
    % get exp decay tau based on different amplitude regions in wave
    [~,maxTprof]= max(deskewedWaveProf);

    % take ~10 ms wide stripe, odd number of pixels
    nPxT = round(10/pxSzT);
    if mod(nPxT,2)~=0
        nPxT = nPxT+1;
    end
    waveProfX = mean(waveImgDeskewed(:,maxTprof-nPxT/2:maxTprof+nPxT/2),2);
  
    % fit with spline
    nK = ceil( (max(x)-min(x))/1 ); % fit every 1 um
    splOrd = 3; % spline order, poly3
    spl0  = spap2(nK, splOrd, x, waveProfX);  % first guess of spline coefficients
    % newknt for a possibly better knot distribution
    knots = newknt(spl0);
    spl0  = spap2(knots, splOrd, x, waveProfX);
    yf_x = fnval(spmak(spl0.knots,spl0.coefs),x);
    
    % find separation between regions with different amplitudes
    [~,minLocs] = findpeaks(-yf_x);
    
    % calculate exp decay time constant for different regions 
    % exp decay function 
    funExpFit = @(p,x,y) (p(1)*exp(-(x-p(2))/p(3))+p(4)) - y;
    funExp = @(p,x) p(1)*exp(-(x-p(2))/p(3))+p(4);
    optFit = optimoptions('lsqnonlin','Display','off');
    
    tauD_regionWise = num2cell(zeros(numel(minLocs)+1,3));
    nRegions = numel(minLocs)+1;
    for i=1:nRegions
        if i==1
            tProfReg = mean(waveImgDeskewed(1:minLocs(i),:),1);
            regPx = [1 minLocs(i)];% region pixels
            regU =  regPx.*pxSzX;% region units, um
        elseif i==nRegions
            tProfReg = mean(waveImgDeskewed(minLocs(i-1):end,:),1);
            regPx = [minLocs(i-1) size(waveImgDeskewed,1)];% region pixels
            regU =  regPx.*pxSzX;% region units, um
        else
            tProfReg = mean(waveImgDeskewed(minLocs(i-1):minLocs(i),:),1);
            regPx = [minLocs(i-1) minLocs(i)];% region pixels
            regU =  regPx.*pxSzX;% region units, um
        end
        
        % fit with spline
        nK = ceil( (max(t)-min(t))/20 ); % fit every 20 ms
        splOrd = 3; % spline order, poly3
        spl0  = spap2(nK, splOrd, t, tProfReg);  % first guess of spline coefficients
        % newknt for a possibly better knot distribution
        knots = newknt(spl0);
        spl0  = spap2(knots, splOrd, t, tProfReg);
        yf = fnval(spmak(spl0.knots,spl0.coefs),t);
        
        % find tr*max in decay part
        [vMax,indMax] = max(yf);
        
        tr = 0.95;
        [~,expDstartFit] = min( abs( yf(indMax:end) - tr*vMax ) );
    
        % exp decay fit 
        t_exp = t(indMax+expDstartFit-1:end);
        y_exp = tProfReg(indMax+expDstartFit-1:end);
        coefExp = lsqnonlin(@(p)funExpFit(p,t_exp,y_exp),...
            [vMax,t_exp(1),20,0],[0 0 0 -inf],[inf inf inf inf],optFit);
        % tauD in ms; position of region in pixels and in um
        tauD_regionWise(i,1) = {coefExp(3)};
        tauD_regionWise(i,2) = {regPx};
        tauD_regionWise(i,3) = {regU};
         
%         figure
%         plot(t,tProfReg,'ok')
%         hold on
%         plot(t,yf,'r')    
%         plot(t_exp,y_exp,'b')
%         plot(t_exp,funExp(coefExp,t_exp),'m','LineWidth',2)
       
    end
    
    
    % clear axes
    delete(hObjsA.ax_detEvent.Children)
    delete(hObjsA.ax_deskewed.Children)
    
    % show analysis
    image(waveImg,'YData',[min(x) max(x)],...
        'XData',[min(t) max(t)],...
        'CDataMapping','scaled','Parent',hObjsA.ax_detEvent);
    set(hObjsA.ax_detEvent,'YGrid','off','FontSize',14)
    set(get(hObjsA.ax_detEvent,'Ylabel'),'String','x (\mum)','FontWeight','bold')
    set(get(hObjsA.ax_detEvent,'Xlabel'),'String','t (ms)','FontWeight','bold')
    set(get(hObjsA.ax_detEvent,'title'),'String','detected event','FontWeight','bold')
    
    image(waveImgDeskewed,'YData',[min(x) max(x)],...
        'XData',[min(t) max(t)],...
        'CDataMapping','scaled','Parent',hObjsA.ax_deskewed);
    set(hObjsA.ax_deskewed,'XTick',[],'YGrid','off','FontSize',14)
    set(get(hObjsA.ax_deskewed,'Ylabel'),'String','x (\mum)','FontWeight','bold')
    set(get(hObjsA.ax_deskewed,'title'),'String','deskewed event','FontWeight','bold')
    
    % set x axis
    hObjsA.ax_detEvent.XLim = [min(t) max(t)];

%     % show different amplitudes regions 
%     tauD_regionWise
%     
%     regLinesY = unique(reshape(regLinesY',[],1));
%   
%     line(hObjsA.ax_deskewed.XLim, [regLinesY regLinesY],...
%         'Parent',hObjsA.ax_deskewed,'Color','k','LineStyle',...
%         '-','LineWidth',1,'Tag','amplRegions');
%     
%     line(yf_x,x,...
%         'Parent',hObjsA.ax_deskewed,'Color','r','LineStyle',...
%         '-','LineWidth',1,'Tag','amplRegions');
%   
%     delete( findall(0,'Type','Line','Tag','amplRegions')  )
    
    % show deskewed profile
    line(t,deskewedWaveProf,...
        'Parent',hObjsA.ax_deskewedProf,'Color',[0 0 0 0.75],'LineStyle','-','LineWidth',1,...
        'Tag','profile');
    
    % save data
    caffeineAnalysis = struct('caffeineImg',waveImg,'caffeineImgDeskewed',waveImgDeskewed,...
                          'tauD_regionWise',{tauD_regionWise},...
                          'deskewedCaffeineProf',deskewedWaveProf);
    
    selectedEvent.analysis.waveAnalysis = [];
    selectedEvent.analysis.transientAnalysis = [];
    selectedEvent.analysis.caffeineAnalysis = caffeineAnalysis;
    selectedEvent.analysis.accepted = false;
    
    setappdata(analysisFig,'selectedEvent',selectedEvent)
    
    % fit t-profile of deskewed wave with spline and with exponential rise
    % to get t0
    findAndAnalyzePeaks([],[],analysisFig)
    
    % set up zoom
    hObjsA.ax_deskewedProf.ButtonDownFcn = @showBiggerAxes;
    
    
   
case 'transient'
    %% %%%%%%%%%%%%%%%%%%%% TRANSIENTS ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % set up zoom
    hObjsA.ax_orgImgProf.ButtonDownFcn = @showBiggerAxes;
    
    % save data
    transientsProf = mean(img,1);
    transientAnalysis = struct('transientsImg',img,...
        'transientsProf',transientsProf);
    
    selectedEvent.analysis.waveAnalysis = [];
    selectedEvent.analysis.transientAnalysis = transientAnalysis;
    selectedEvent.analysis.caffeineAnalysis = [];
    selectedEvent.analysis.accepted = false;
    
    setappdata(analysisFig,'selectedEvent',selectedEvent)
    
    % fit t-profile of image with spline and with exponential rise
    % to get t0
    findAndAnalyzePeaks([],[],analysisFig)
    
    
end

set(analysisFig,'Pointer','arrow')
drawnow
end










