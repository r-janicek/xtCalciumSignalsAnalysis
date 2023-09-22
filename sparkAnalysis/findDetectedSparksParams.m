function eventParams = findDetectedSparksParams(...
    img, statSparks, mainFig ,calcMethod, sparksRec)

% img = not normalized data, raw data
% statSparks = statistic from sparks detection

% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

pxSzT = imgData.pxSzT;
pxSzX = imgData.pxSzX;
scrRes = get(0,'ScreenSize');

% filter raw image
img = imgFiltering(img, pxSzT, pxSzX);

% t-profile = average from 1 um width area
% x-profile = average from 5 ms width area
n_px_t = ceil(1/pxSzX);
if mod(n_px_t,2)==0
    n_px_t = n_px_t - 1;
end

n_px_x = ceil(5/pxSzT);
if mod(n_px_x,2)==0
    n_px_x = n_px_x - 1;
end

% out = eventBoundaries(img,statEvents,mainFig)
% change pointer to watch
set(mainFig,'Pointer','watch')
drawnow

% check if there was some sparks detected
if ~isempty(statSparks)
    % add waitbar
    hw = waitbar(0,'0 %','Name','analyzing detected sparks...',...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(hw,'canceling',0)
    
    if isfield(getappdata(mainFig),'batchProcessing')
        if getappdata(mainFig,'batchProcessing')
            hw.Visible = 'off';
        end
    end
    
    for i = 1:length(statSparks)

        % get number of spark ROI rectangle
        sparkROInum = str2double(sparksRec(i).Tag);
        
        % get image of event and max crossing profiles
        [imgE, imgEs, imgE_m, maxOfEventPos, ... 
            maxCrossProfs, t0, bs, eventROIstart_t] = ...
            eventImgAndMaxCrossingProfiles( ...
            statSparks(i), [], [],...
            pxSzT, pxSzX, n_px_t, n_px_x, img, [], mainFig);

        t = maxCrossProfs.t;
        t_ups = maxCrossProfs.t_ups;
        t_event_prof = maxCrossProfs.t_event_prof;
        t_event_profS = maxCrossProfs.t_event_prof_S;
        
        peakData_tProf = maxCrossProfs.peakData_tProf;
        
        x = maxCrossProfs.x;
        x_ups = maxCrossProfs.x_ups;
        x_event_prof = maxCrossProfs.x_event_prof;
        x_event_profS = maxCrossProfs.x_event_prof_S;
        
        r_m = maxOfEventPos(1);
        c_m = maxOfEventPos(2);
        % masks of profiles
        tProf_m = maxCrossProfs.t_event_prof_m;
        xProf_m = imgE_m(:,c_m);
            
        % choose method how to calculate spark parameters
        calcMethodPrev = calcMethod;
        switch calcMethod
            
            case '2DGauss'
                % calculate parameters, fit spark data with 2D exponentially modified gaussian
                try
                    try
                        % create time and spatial grids
                        [T,X] = meshgrid(t,x);
                        [T_ups,X_ups] = meshgrid(t_ups,x_ups);
                        
                        gaussModel = '2D_EMGxt';
                        fit2DGauss = fit2DGaussian(...
                            T,X,T_ups,X_ups,imgE,imgE_m,maxOfEventPos,gaussModel);
                        
                        imgE_fit = fit2DGauss.dataEventFit;
                        imgE_fit_ups = fit2DGauss.dataEventFit_ups;
                        
                        % get t and x profiles crossing max of fit
                        [rFitMax,cFitMax] = find(imgE_fit == max(imgE_fit(:)));
                        t_data_prof = imgE(rFitMax,:); t_data_prof = t_data_prof(:);
                        x_data_prof = imgE(:,cFitMax); x_data_prof = x_data_prof(:);
                        [rFitMaxUps,cFitMaxUps] = find(imgE_fit_ups == max(imgE_fit_ups(:)));
                        t_fit_prof = imgE_fit_ups(rFitMaxUps,:);
                        x_fit_prof = imgE_fit_ups(:,cFitMaxUps);
                        t_fit_prof = t_fit_prof(:);
                        x_fit_prof = x_fit_prof(:);
                        
                        [~,peakData_ups.pos] = ...
                            min( abs( t_ups - t(maxCrossProfs.peakData_tProf.pos) ) );
                        peakData_ups.val = t_fit_prof(peakData_ups.pos);
                        
                        eP_2Dfit = getParametersOfEventProfile(t_ups,t_fit_prof,...
                            x_ups,x_fit_prof,...
                            t_fit_prof(1),t_fit_prof(1),[],...
                            imgData.blank,peakData_ups,[]);
                        
                    catch
                        % get params from image of event
                        calcMethodPrev = '2DGauss'; 
                        calcMethod = 'estimate from event img';
                        
                        eP_2Dfit = getParametersOfEventFromRawImage(...
                            imgE, imgEs, imgE_m, ...
                            [], [], r_m, c_m, ...
                            t_event_profS, x_event_profS, t, x, imgData);

                    end
                    
                    eventParams.amplitude(i,1) = eP_2Dfit.amplitude;
                    eventParams.t0(i,1) = eP_2Dfit.t0 + eventROIstart_t;
                    eventParams.TTP(i,1) = eP_2Dfit.TTP;
                    eventParams.FDHM(i,1) = eP_2Dfit.FDHM;
                    eventParams.FWHM(i,1) = eP_2Dfit.FWHM;
                    eventParams.sparkMass(i,1) = eP_2Dfit.sparkMass;
                    
                    try
                        % only if fit by EMG
                        eventParams.tauD(i,1) = fit2DGauss.coef2DGauss(4);
                        % (F-F0)/(F0-blank)
                        imgE_fit_ups_norm = ...
                            (imgE_fit_ups - fit2DGauss.coef2DGauss(end))/(fit2DGauss.coef2DGauss(end) - imgData.blank);
                        % area under curve
                        eventParams.AUC_2DFit(i,1) = ...
                            sum(imgE_fit_ups_norm(:).* mean(diff(t_ups)).*mean(diff(x_ups)));
                        %eP_2Dfit.AUC_2DFit = trapz(x_ups,trapz(t_ups, D_fit_ups_norm, 2));
                    catch
                        eventParams.AUC_2DFit(i,1) = nan;
                    end
                     
                catch
                    eventParams.amplitude(i,1) = nan;
                    eventParams.t0(i,1) = nan;
                    eventParams.TTP(i,1) = nan;
                    eventParams.FDHM(i,1) = nan;
                    eventParams.FWHM(i,1) = nan;
                    eventParams.sparkMass(i,1) = nan;
                    eventParams.tauD(i,1) = nan;
                    eventParams.AUC_2DFit(i,1) = nan;
                end
                  
            case 'peakXTProfiles'
                      
                try
                    % get parameters of spark from profiles without 2D fitting            
                    % do fits of profiles
                    try
                        % try fitting
                        % '1expR1expD' piecewise: coefficients: [t0, F01, tauR, A, t1, tauD, F02]
                        % 'spline' piecewise: coefficients: [t0, F01]
                        % modified gaussian = 'EMG', coefficients:[A,m,sd,tau,F0]
                        % p0 = [t0 bs tauR v_max_t p_max_t*pxSzT 15 mean(t_spark_prof(end-5:end))];
                        % p0 = [v_max_t p_max_t*pxSzT p_max_t*pxSzT-t0 15 bs];
                        % initial parameters for time profile fit
                        p0 = [t0 bs];
                        out_t_prof = fitOneWholeEvent( ...
                            p0, t, t_ups, t_event_prof, ...
                            'yes', 'spline', tProf_m);
                        
                        bs_t_prof = out_t_prof.bs;
                        t0_t_prof = out_t_prof.t0;
                        tauD = out_t_prof.tauD;
                        
                        % fit gaussian or exp modified gaussian or spline
                        % to spatial profile
                        % 'Gauss' : [F0,A,w,xc]
                        % 'EMG' : [A,m,sd,tau,F0]
                        x_event_prof_m = x_event_prof;
                        x_event_prof_m(~xProf_m) = nan;
                        [v_max_x,p_max_x] = max(x_event_prof_m);
                        
                        % estimation of the width in half max
                        w = sum((x_event_prof_m-bs-(v_max_x-bs)/2)>=0)*pxSzX;
                        % fit
                        out_x_prof = fitOneWholeEvent(...
                            [bs v_max_x-bs w p_max_x*pxSzX-pxSzX],...
                            x, x_ups, x_event_prof, 'no', 'Gauss', xProf_m);
                        
                        % baseline of x profile
                        bs_x_prof = out_x_prof.bs;
                        
                        % upscale mask
                        % start and end of event's mask
                        tProf_m_ups = false(size(t_ups));
                        sE = find(tProf_m,1,'first');
                        eE = find(tProf_m,1,'last');
                        [~, sE_ups] = min( abs( t_ups - t(sE) ) );
                        [~, eE_ups] = min( abs( t_ups - t(eE) ) );
                        tProf_m_ups(sE_ups:eE_ups) = true;
                        % peak position and value
                        [~, loc_ups] = findpeaks(...
                            out_t_prof.yFit_ups(tProf_m_ups), ...
                            t_ups(tProf_m_ups), 'NPeaks',1, ...
                            'SortStr','descend');
                        [~, peakData_ups.pos] = ...
                            min( abs( t_ups - loc_ups ) );
                        peakData_ups.val = ...
                            out_t_prof.yFit_ups(peakData_ups.pos);

                        eP_profs = getParametersOfEventProfile( ...
                            t_ups,out_t_prof.yFit_ups,...
                            x_ups,out_x_prof.yFit_ups,...
                            bs_t_prof, bs_x_prof, t0_t_prof,...
                            imgData.blank, peakData_ups, tProf_m_ups);
                       
                    catch
                        % get params from image of event
                        calcMethodPrev = 'peakXTProfiles'; 
                        calcMethod = 'estimate from event img';
                        
                        eP_profs = getParametersOfEventFromRawImage(...
                            imgE, imgEs, imgE_m, ...
                            [], [], r_m, c_m, ...
                            t_event_profS, x_event_profS, t, x, imgData);
                    end
                    
                    eventParams.amplitude(i,1) = eP_profs.amplitude;
                    eventParams.t0(i,1) = eP_profs.t0 + eventROIstart_t;
                    eventParams.TTP(i,1) = eP_profs.TTP;
                    eventParams.FDHM(i,1) = eP_profs.FDHM;
                    eventParams.FWHM(i,1) = eP_profs.FWHM;
                    eventParams.sparkMass(i,1) = eP_profs.sparkMass;
                    eventParams.tauD(i,1) = tauD;
                    eventParams.AUC_2DFit(i,1) = nan;
                      
                catch
                 
                    % if cannot get profiles for some reason
                    eventParams.amplitude(i,1) = nan;
                    eventParams.t0(i,1) = nan;
                    eventParams.TTP(i,1) = nan;
                    eventParams.FDHM(i,1) = nan;
                    eventParams.FWHM(i,1) = nan;
                    eventParams.sparkMass(i,1) = nan;
                    eventParams.tauD(i,1) = nan;
                    eventParams.AUC_2DFit(i,1) = nan;
                    
                end
        
                
            case 'estimate from event img'
                try
                    % get parameters of spark from profiles without 2D fitting
                    % use coeficient from previous fitting (t0, tauR, A, baseline)
                    % get params from image of event
                    if ~isempty(t0)
                        eParamsEst = ...
                            getParametersOfEventFromRawImage(...
                            imgE, imgEs, imgE_m, ...
                            bs, t0, r_m, c_m, ...
                            t_event_profS, x_event_profS, t, x, imgData);
                    else
                        eParamsEst = ...
                            getParametersOfEventFromRawImage(...
                            imgE, imgEs, imgE_m, ...
                            [], [], r_m, c_m, ...
                            t_event_profS, x_event_profS, t, x, imgData);
                    end
                    
                    eventParams.amplitude(i,1) = eParamsEst.amplitude;
                    eventParams.t0(i,1) = eParamsEst.t0 + eventROIstart_t;
                    eventParams.TTP(i,1) = eParamsEst.TTP;
                    eventParams.FDHM(i,1) = eParamsEst.FDHM;
                    eventParams.FWHM(i,1) = eParamsEst.FWHM;
                    eventParams.sparkMass(i,1) = eParamsEst.sparkMass;
                    eventParams.tauD(i,1) = eParamsEst.tauD;
                    eventParams.AUC_2DFit(i,1) = nan;
                    
                catch
                    % if cannot get estimate of params for some reason
                    eventParams.amplitude(i,1) = nan;
                    eventParams.t0(i,1) = nan;
                    eventParams.TTP(i,1) = nan;
                    eventParams.FDHM(i,1) = nan;
                    eventParams.FWHM(i,1) = nan;
                    eventParams.sparkMass(i,1) = nan;
                    eventParams.tauD(i,1) = nan;
                    eventParams.AUC_2DFit(i,1) = nan;
                end
                
        end
        
        eventParams.sparkROINum(i,1) = sparkROInum;
        
        % show  individual calcium events
        if hObjs.check_showEventsFigs.Value
            % create figure
            scF = 0.8;
            hf = figure('Tag','CaEventParam',...
                'Name',sprintf('spark from ROI #: %d',sparkROInum),...
                'Visible','on',...
                'OuterPosition',[scrRes(3)-scrRes(3)*scF ...
                scrRes(4)*0.4 ...
                scrRes(3)*scF ...
                scrRes(4)*0.6]);
            set(hf,'Units','Normalized')
            
            h_a = 0.55;
            w_a = 0.225;
            dx = 0.05;
            dy = 0.07;
            
            fontSzT = 18;
            fontSzL = 16;
            fontSzNum = 14;
            fontSzLegend = 14;
            
            switch calcMethod
                
                case '2DGauss'
                    try
                        %%%%% 2D gauss fit params %%%%%
                        uicontrol('Style','text','Parent',hf,'Units','normalized',...
                            'Position', [dx 1-dy 1-2*dx dy],'FontUnits','normalized','FontSize',0.4,...
                            'FontWeight','bold','HorizontalAlignment','left',...
                            'String',sprintf('parameters from 2D exponentially modified gaussian fit (%s)',gaussModel));
                        
                        %
                        ha1 = axes(hf,'Units','normalized','Position',[dx 1-h_a-2*dy w_a*1.5 h_a]);
                        mesh(T,X,imgE,'LineStyle','-','LineWidth',0.5,'FaceColor','none','EdgeColor','k','EdgeAlpha',0.4,'Parent',ha1)
                        %plot3(T,X,D,'LineStyle','-','Marker','.','Color',[0.1 0.1 0.1],'Parent',ha1)
                        hold on
                        surf(T,X,imgE_fit,'FaceAlpha',0.6,'EdgeColor','none','FaceColor','interp','Parent',ha1)
                        line(t_ups,X_ups(rFitMaxUps,:),t_fit_prof,'Parent',ha1,'LineWidth',3,'Color','k')
                        line(T_ups(:,cFitMaxUps),x_ups,x_fit_prof,'Parent',ha1,'LineWidth',3,'Color','k')
                        set(ha1,'XLim',[min(t) max(t)],'YLim',[min(x) max(x)],'ZLim',[min([imgE(:);imgE_fit(:)])*0.95 max([imgE(:);imgE_fit(:)])*1.05])
                        colormap(jet)
                        view(-15,40)
                        set(ha1,'YDir','normal','FontSize',fontSzNum)
                        title(ha1,sprintf('spark from ROI #: %d -- 2D gauss fit',sparkROInum),'FontSize',fontSzT)
                        xlabel(ha1,'t (ms)','FontSize',fontSzL)
                        ylabel(ha1,'x (\mum)','FontSize',fontSzL)
                        zlabel(hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
                  
                        %
                        ha2 = axes(hf,'Units','normalized','Position',[ha1.Position(1)+ha1.Position(3)+dx 1-h_a-2*dy w_a h_a]);
                        line(x,x_data_prof,'Parent',ha2,'LineWidth',1,'Color','k','Display','data')
                        line(x_ups,x_fit_prof,'Parent',ha2,'LineWidth',2,'Color','r','Display','fit')
                        set(ha2,'FontSize',fontSzNum)
                        title(ha2,sprintf('spark from ROI #: %d -- x profile',sparkROInum),'FontSize',fontSzT)
                        ylabel(ha2,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
                        xlabel(ha2,'x (\mum)','FontSize',fontSzL)
                        xlim(ha2,[min(x_ups) max(x_ups)])
                        ylim(ha2,[min([x_data_prof;x_fit_prof])*0.95 max([x_data_prof;x_fit_prof])*1.05])
                        line([eP_2Dfit.half_max_x_1 eP_2Dfit.half_max_x_2],[eP_2Dfit.half_max_x eP_2Dfit.half_max_x],'Parent',ha2,...
                            'LineWidth',2,'Color','b','Display','FWHM')
                        hl_2 = legend(ha2,'show');
                        hl_2.Location = 'best';
                        hl_2.FontSize = fontSzLegend;
                        
                        %
                        ha3 = axes(hf,'Units','normalized','Position',[ha2.Position(1)+ha2.Position(3)+dx 1-h_a-2*dy w_a h_a]);
                        line(t,t_data_prof,'Parent',ha3,'LineWidth',1,'Color','k','Display','data')
                        line(t_ups,t_fit_prof,'Parent',ha3,'LineWidth',2,'Color','r','Display','fit')
                        set(ha3,'FontSize',fontSzNum)
                        title(ha3,sprintf('spark from ROI #: %d -- t profile',sparkROInum),'FontSize',fontSzT)
                        xlabel(ha3,'t (ms)','FontSize',fontSzL)
                        ylabel(ha3,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
                        xlim(ha3,[min(t_ups) max(t_ups)])
                        ylim(ha3,[min([t_data_prof;t_fit_prof])*0.95 max([t_data_prof;t_fit_prof])*1.05])
                        line([eP_2Dfit.half_max_t_1 eP_2Dfit.half_max_t_2],[eP_2Dfit.half_max_t eP_2Dfit.half_max_t],'Parent',ha3,...
                            'LineWidth',2,'Color','b','Display','FDHM')
                        line([eP_2Dfit.t_max eP_2Dfit.t_max],[eP_2Dfit.bs_t eP_2Dfit.v_max],'Parent',ha3,'LineWidth',2,'Color','g','Display','amplitude')
                        line([eP_2Dfit.t0 eP_2Dfit.t_max],[eP_2Dfit.bs_t eP_2Dfit.bs_t],'Parent',ha3,'LineWidth',2,'Color','m','Display','TTP')
                        hl_3 = legend(ha3,'show');
                        hl_3.Location = 'best';
                        hl_3.FontSize = fontSzLegend;
                        
                        %
                        uicontrol('Style','text','Parent',hf,'Units','normalized',...
                            'Position', [dx 1-4*dy-h_a-2*dy 1-2*dx 2*dy],'FontUnits','normalized','FontSize',0.25,...
                            'FontWeight','normal','HorizontalAlignment','left',...
                            'String',[...
                            sprintf('parameters of spark from ROI #: (#%d):',sparkROInum),' // ',...
                            sprintf('amplitude = %0.2f ',eP_2Dfit.amplitude),'(',char(916),'F/F0)',' // ',...
                            sprintf('TTP = %0.2f (ms)',eP_2Dfit.TTP),' // ',...
                            sprintf('FDHM = %0.2f (ms)',eP_2Dfit.FDHM),' // ',...
                            sprintf('FWHM = %0.2f ',eP_2Dfit.FWHM),'(',char(181),'m)',' // ',...
                            sprintf('sparkMass = %0.2f ',eP_2Dfit.sparkMass),'(',char(916),'F/F0*',char(181),'m^3)' ]);
                    catch
                        
                    end
                    
                case 'peakXTProfiles'
                    try
                        %%%%% profiles fit params %%%%%
                        uicontrol('Style','text','Parent',hf,'Units','normalized',...
                            'Position', [dx 1-dy 1-2*dx dy],'FontUnits','normalized','FontSize',0.4,...
                            'FontWeight','bold','HorizontalAlignment','left',...
                            'String','parameters from profiles (1 um and 5 ms wide) crossing peak of event');
                        %
                        ha1 = axes(hf,'Units','normalized','Position',[dx 1-h_a-2*dy w_a*1.5 h_a]);
                        image(imgE,'YData',[min(x_ups) max(x_ups)],'XData',[min(t_ups) max(t_ups)],'CDataMapping','scaled','Parent',ha1)
                        line(get(gca,'XLim'),[(r_m-1-(n_px_t-1)/2)*pxSzX (r_m-1-(n_px_t-1)/2)*pxSzX],'Parent',ha1,'LineWidth',1,'Color','k',...
                            'LineStyle','-');
                        line(get(gca,'XLim'),[(r_m-1+(n_px_t-1)/2)*pxSzX (r_m-1+(n_px_t-1)/2)*pxSzX],'Parent',ha1,'LineWidth',1,'Color','k',...
                            'LineStyle','-');
                        
                        line([(c_m-1-(n_px_x-1)/2)*pxSzT (c_m-1-(n_px_x-1)/2)*pxSzT],get(gca,'YLim'),'Parent',ha1,'LineWidth',1,'Color','k',...
                            'LineStyle','-');
                        line([(c_m-1+(n_px_x-1)/2)*pxSzT (c_m-1+(n_px_x-1)/2)*pxSzT],get(gca,'YLim'),'Parent',ha1,'LineWidth',1,'Color','k',...
                            'LineStyle','-');
                        set(ha1,'FontSize',fontSzNum)
                        title(ha1,sprintf('spark from ROI #: %d',sparkROInum),'FontSize',fontSzT)
                        xlabel(ha1,'t (ms)','FontSize',fontSzL)
                        ylabel(ha1,'x (\mum)','FontSize',fontSzL)
                        
                        %
                        ha2 = axes(hf,'Units','normalized','Position',[ha1.Position(1)+ha1.Position(3)+dx 1-h_a-2*dy w_a h_a]);
                        line(x,x_event_prof,'Parent',ha2,'LineWidth',1,'Color','k','Display','data')
                        line(x_ups,out_x_prof.yFit_ups,'Parent',ha2,'LineWidth',2,'Color','r',...
                            'Display',sprintf('fit with %s',out_x_prof.fitModel))
                        set(ha2,'FontSize',fontSzNum)
                        title(sprintf('spark from ROI #: %d -- x profile',sparkROInum),'FontSize',fontSzT)
                        ylabel(ha2,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
                        xlabel(ha2,'x (\mum)','FontSize',fontSzL)
                        xlim(ha2,[min(x_ups) max(x_ups)])
                        ylim(ha2,[min([x_event_prof;out_x_prof.yFit_ups])*0.95 max([x_event_prof;out_x_prof.yFit_ups])*1.05])
                        line([eP_profs.half_max_x_1 eP_profs.half_max_x_2],[eP_profs.half_max_x eP_profs.half_max_x],...
                            'Parent',ha2,'LineWidth',2,'Color','b','Display','FWHM')
                        hl_2 = legend(ha2,'show');
                        hl_2.Location = 'best';
                        hl_2.FontSize = fontSzLegend;
                        
                        %
                        ha3 = axes(hf,'Units','normalized','Position',[ha2.Position(1)+ha2.Position(3)+dx 1-h_a-2*dy w_a h_a]);
                        line(t,t_event_prof,'Parent',ha3,'LineWidth',1,'Color','k','Display','data')
                        line(t_ups,out_t_prof.yFit_ups,'Parent',ha3,'LineWidth',2,'Color','r',...
                            'Display',sprintf('fit with %s',out_t_prof.fitModel))
                        set(ha3,'FontSize',fontSzNum)
                        title(ha3,sprintf('spark from ROI #: %d -- t profile',sparkROInum),'FontSize',fontSzT)
                        xlabel(ha3,'t (ms)','FontSize',fontSzL)
                        ylabel(ha3,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
                        xlim(ha3,[min(t_ups) max(t_ups)])
                        ylim(ha3,[min([t_event_prof;out_t_prof.yFit_ups])*0.95 max([t_event_prof;out_t_prof.yFit_ups])*1.05])
                        line([eP_profs.half_max_t_1 eP_profs.half_max_t_2],[eP_profs.half_max_t eP_profs.half_max_t],...
                            'Parent',ha3,'LineWidth',2,'Color','b','Display','FDHM')
                        line([eP_profs.t_max eP_profs.t_max],[eP_profs.bs_t eP_profs.v_max],'Parent',ha3,'LineWidth',2,'Color','g','Display','amplitude')
                        line([eP_profs.t0 eP_profs.t_max],[eP_profs.bs_t eP_profs.bs_t],'Parent',ha3,'LineWidth',2,'Color','m','Display','TTP')
                        hl_3 = legend(ha3,'show');
                        hl_3.Location = 'best';
                        hl_3.FontSize = fontSzLegend;
                        
                        %
                    
                        uicontrol('Style','text','Parent',hf,'Units','normalized',...
                            'Position', [dx 1-4*dy-h_a-2*dy 1-2*dx 2*dy],'FontUnits','normalized','FontSize',0.25,...
                            'FontWeight','normal','HorizontalAlignment','left',...
                            'String',[...
                            sprintf('parameters of spark from ROI #: (#%d):',sparkROInum),' // ',...
                            sprintf('amplitude = %0.2f ',eP_profs.amplitude),'(',char(916),'F/F0)',' // ',...
                            sprintf('TTP = %0.2f (ms)',eP_profs.TTP),' // ',...
                            sprintf('FDHM = %0.2f (ms)',eP_profs.FDHM),' // ',...
                            sprintf('FWHM = %0.2f ',eP_profs.FWHM),'(',char(181),'m)',' // ',...
                            sprintf('sparkMass = %0.2f ',eP_profs.sparkMass),'(',char(916),'F/F0*',char(181),'m^3)' ]);
                        
                    catch
                        
                    end
                   
                case 'estimate from event img'
                    try
                        if ~strcmp(calcMethod,calcMethodPrev)
                            eParamsEst = eP_profs;
                        end

                        %%%%% from raw image params %%%%%
                        uicontrol('Style','text','Parent',hf,'Units','normalized',...
                            'Position', [dx 1-dy 1-2*dx dy],'FontUnits','normalized','FontSize',0.4,...
                            'FontWeight','bold','HorizontalAlignment','left',...
                            'String', 'parameters from profiles (1 um and 5 ms wide) crossing peak of event (from image smoothed with 2D spline )');
                        %
                        ha1 = axes(hf,'Units','normalized','Position',[dx 1-h_a-2*dy w_a*1.5 h_a]);
                        image(imgE,'YData',[min(x_ups) max(x_ups)],'XData',[min(t_ups) max(t_ups)],'CDataMapping','scaled','Parent',ha1)
                        line(get(gca,'XLim'),[(r_m-1-(n_px_t-1)/2)*pxSzX (r_m-1-(n_px_t-1)/2)*pxSzX],'Parent',ha1,'LineWidth',1,'Color','k',...
                            'LineStyle','-');
                        line(get(gca,'XLim'),[(r_m-1+(n_px_t-1)/2)*pxSzX (r_m-1+(n_px_t-1)/2)*pxSzX],'Parent',ha1,'LineWidth',1,'Color','k',...
                            'LineStyle','-');
                        
                        line([(c_m-1-(n_px_x-1)/2)*pxSzT (c_m-1-(n_px_x-1)/2)*pxSzT],get(gca,'YLim'),'Parent',ha1,'LineWidth',1,'Color','k',...
                            'LineStyle','-');
                        line([(c_m-1+(n_px_x-1)/2)*pxSzT (c_m-1+(n_px_x-1)/2)*pxSzT],get(gca,'YLim'),'Parent',ha1,'LineWidth',1,'Color','k',...
                            'LineStyle','-');
                        set(ha1,'FontSize',fontSzNum)
                        title(ha1,sprintf('spark from ROI #: %d',sparkROInum),'FontSize',fontSzT)
                        xlabel(ha1,'t (ms)','FontSize',fontSzL)
                        ylabel(ha1,'x (\mum)','FontSize',fontSzL)
                 
                        %
                        ha2 = axes(hf,'Units','normalized','Position',[ha1.Position(1)+ha1.Position(3)+dx 1-h_a-2*dy w_a h_a]);
                        line(x,x_event_prof,'Parent',ha2,'LineWidth',1,'Color','k','Display','data')
                        line(x,x_event_profS,'Parent',ha2,'LineWidth',2,'Color','r',...
                            'Display','profile from image smoothed with 2D spline')
                        set(ha2,'FontSize',fontSzNum)
                        title(sprintf('spark from ROI #: %d -- x profile',sparkROInum),'FontSize',fontSzT)
                        ylabel(ha2,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
                        xlabel(ha2,'x (\mum)','FontSize',fontSzL)
                        xlim(ha2,[min(x) max(x)])
                        ylim(ha2,[ min([x_event_prof;x_event_profS])*0.95 ...
                            max([x_event_prof;x_event_profS])*1.05])
                        line([eParamsEst.half_max_x_1 eParamsEst.half_max_x_2], ...
                            [eParamsEst.half_max_x eParamsEst.half_max_x],...
                            'Parent',ha2,'LineWidth',2,'Color','b','Display','FWHM')
                        hl_2 = legend(ha2,'show');
                        hl_2.Location = 'best';
                        hl_2.FontSize = fontSzLegend;
                        
                        %
                        ha3 = axes(hf,'Units','normalized','Position',[ha2.Position(1)+ha2.Position(3)+dx 1-h_a-2*dy w_a h_a]);
                        line(t,t_event_prof,'Parent',ha3,'LineWidth',1,'Color','k','Display','data')
                        line(t,t_event_profS,'Parent',ha3,'LineWidth',2,'Color','r',...
                            'Display','profile from image smoothed with 2D spline')
                        set(ha3,'FontSize',fontSzNum)
                        title(ha3,sprintf('spark from ROI #: %d -- t profile',sparkROInum),'FontSize',fontSzT)
                        xlabel(ha3,'t (ms)','FontSize',fontSzL)
                        ylabel(ha3,hObjs.ax_prof.YLabel.String,'FontSize',fontSzL)
                        xlim(ha3,[min(t) max(t)])
                        ylim(ha3,[min([t_event_prof;t_event_profS])*0.95 max([t_event_prof;t_event_profS])*1.05])
                        line([eParamsEst.half_max_t_1 eParamsEst.half_max_t_2],...
                            [eParamsEst.half_max_t eParamsEst.half_max_t],...
                            'Parent',ha3,'LineWidth',2,'Color','b','Display','FDHM')
                        line([eParamsEst.t_max eParamsEst.t_max],...
                            [eParamsEst.bs eParamsEst.v_max],...
                            'Parent',ha3,'LineWidth',2,'Color','g','Display','amplitude')
                        line([eParamsEst.t0 eParamsEst.t_max],[eParamsEst.bs eParamsEst.bs],'Parent',ha3,'LineWidth',2,'Color','m','Display','TTP')
                        hl_3 = legend(ha3,'show');
                        hl_3.Location = 'best';
                        hl_3.FontSize = fontSzLegend;
                        
                        %
                        uicontrol('Style','text','Parent',hf,'Units','normalized',...
                            'Position', [dx 1-4*dy-h_a-2*dy 1-2*dx 2*dy],'FontUnits','normalized','FontSize',0.25,...
                            'FontWeight','normal','HorizontalAlignment','left',...
                            'String',[...
                            sprintf('parameters of spark from ROI #: (#%d):',sparkROInum),' // ',...
                            sprintf('amplitude = %0.2f ',eParamsEst.amplitude),'(',char(916),'F/F0)',' // ',...
                            sprintf('TTP = %0.2f (ms)',eParamsEst.TTP),' // ',...
                            sprintf('FDHM = %0.2f (ms)',eParamsEst.FDHM),' // ',...
                            sprintf('FWHM = %0.2f ',eParamsEst.FWHM),'(',char(181),'m)',' // ',...
                            sprintf('sparkMass = %0.2f ',eParamsEst.sparkMass),'(',char(916),'F/F0*',char(181),'m^3)']);
                        
                    catch
                        
                    end
                    
            end
            
        end
        
        % if method for parameters estimation was changed,
        % due to error in selected method, 
        % then after using simpler method for current iteration calculation
        % change it back to selected one
        calcMethod = calcMethodPrev;
        
        clearvars D t x t_prof x_prof t_ups x_ups t_prof_ups x_prof_ups y_t1 y_t2 D1 y_x1 y_x2
        
        clearvars rows_e cols_e eP_profs eP_2Dfit out_t_prof out_x_prof p0 t_spark_prof x_spark_prof
        
        % Report current estimate in the waitbar's message field
        waitbar(i/length(statSparks),hw,sprintf('%d %%',round(100*i/length(statSparks))))
        
    end
    
    % delete waitbar
    delete(hw), clearvars hw
    
else
    eventParams.amplitude = [];
    eventParams.t0 = [];
    eventParams.TTP = [];
    eventParams.FDHM = [];
    eventParams.FWHM = [];
    eventParams.sparkMass = [];
    eventParams.tauD = [];
    eventParams.AUC_2DFit = [];
    eventParams.sparkROINum = [];
end

% put back pointer to arrow
set(mainFig,'Pointer','arrow')
drawnow


end



