function showSelectedEventInBrowser(h_sld, ~, hf_evntBrowser, mainFig)

% get data
imgData = getappdata(mainFig, 'imgData');
analyzedEvntsBrowserTbl = ...
    getappdata(mainFig, 'sparkDetection').analyzedEvntsBrowserTbl;
hObjs_evntsBrowser = getappdata(hf_evntBrowser, 'hObjs');
hObjs_mainFig = getappdata(mainFig, 'hObjs');
% selected event
evnt = analyzedEvntsBrowserTbl(round(h_sld.Value), :);
% set fonts sizes
fontSzNum = 14;
fontSzLegend = 14;
fontSzT = 16;
fontSzL = 18;
% clear axes
ha_img = hObjs_evntsBrowser.ha_img;
ha_xProf = hObjs_evntsBrowser.ha_xProf;
ha_tProf = hObjs_evntsBrowser.ha_tProf;
cla([ha_img, ha_xProf, ha_tProf])
% update slider string
hObjs_evntsBrowser.h_txt_evntSld.String = evnt.axDesc{1}.fig_title_txt;

% show selected event
% show image of event
switch evnt.calcMethod{1}
    case '2DGauss'
        mesh(evnt.axData{1}.T, evnt.axData{1}.X, evnt.axData{1}.imgE, ...
            'LineStyle','-', 'LineWidth',0.5, ...
            'FaceColor','none', 'EdgeColor','k', ...
            'EdgeAlpha',0.4, 'Parent',ha_img)
        %plot3(T, X, D, 'LineStyle','-',
        % 'Marker','.', 'Color',[0.1 0.1 0.1], 'Parent',ha1)
        hold(ha_img, 'on')
        surf(evnt.axData{1}.T, evnt.axData{1}.X, evnt.axData{1}.imgE_fit, ...
            'FaceAlpha',0.6, 'EdgeColor','none', ...
            'FaceColor','interp', 'Parent',ha_img)
        line(evnt.axData{1}.t_ups(:), evnt.axData{1}.x_ups_peak(:), ...
            evnt.axData{1}.t_prof_fit(:), ...
            'Parent',ha_img, 'LineWidth',3, 'Color','k')
        line(evnt.axData{1}.t_ups_peak(:), evnt.axData{1}.x_ups(:), ...
            evnt.axData{1}.x_prof_fit(:), ...
            'Parent',ha_img, 'LineWidth',3, 'Color','k')
        set(ha_img, 'XLim',[min(evnt.axData{1}.t) max(evnt.axData{1}.t)], ...
            'YLim',[min(evnt.axData{1}.x) max(evnt.axData{1}.x)], ...
            'ZLim',[min(evnt.axData{1}.imgE, [],'all'), ...
                    max(evnt.axData{1}.imgE(evnt.axData{1}.imgE_m), [],'all')*1.05])
        old_img_cLims = get(ha_img, 'Clim');
        % in case taht in image of events there is also a part
        % of another event with much higher amplitude
        set(ha_img, 'Clim', [old_img_cLims(1) ...
            prctile(evnt.axData{1}.imgE(evnt.axData{1}.imgE_m),99,'all')])
        colormap(jet)
        view(ha_img, -15,40)
        set(ha_img,'YDir','normal', 'FontSize',fontSzNum)
        title(ha_img, [evnt.axDesc{1}.fig_title_txt,' -- 2D gauss fit'], ...
            'FontSize',fontSzT) 
        xlabel(ha_img,'t (ms)', 'FontSize',fontSzL)
        ylabel(ha_img,'x (\mum)', 'FontSize',fontSzL)
        zlabel(ha_img, hObjs_mainFig.ax_prof.YLabel.String, ...
            'FontSize',fontSzL)
    otherwise
        image(evnt.axData{1}.imgE, ...
            'YData',[min(evnt.axData{1}.x_ups) max(evnt.axData{1}.x_ups)], ...
            'XData',[min(evnt.axData{1}.t_ups) max(evnt.axData{1}.t_ups)], ...
            'CDataMapping','scaled', 'Parent',ha_img)
        % show lines of areas from where profiles are
        % calculated
        r_m = evnt.axData{1}.r_m;
        n_px_t = evnt.axData{1}.n_px_t;
        c_m = evnt.axData{1}.c_m;
        n_px_x = evnt.axData{1}.n_px_x;
        line(get(ha_img,'XLim'), ...
            [(r_m-1-(n_px_t-1)/2)*imgData.pxSzX (r_m-1-(n_px_t-1)/2)*imgData.pxSzX], ...
            'Parent',ha_img, 'LineWidth',1, ...
            'Color','k', 'LineStyle','-');
        line(get(ha_img,'XLim'), ...
            [(r_m-1+(n_px_t-1)/2)*imgData.pxSzX (r_m-1+(n_px_t-1)/2)*imgData.pxSzX], ...
            'Parent',ha_img, 'LineWidth',1, ...
            'Color','k', 'LineStyle','-');
        line([(c_m-1-(n_px_x-1)/2)*imgData.pxSzT (c_m-1-(n_px_x-1)/2)*imgData.pxSzT], ...
            get(ha_img,'YLim'), 'Parent',ha_img, ...
            'LineWidth',1, 'Color','k', 'LineStyle','-');
        line([(c_m-1+(n_px_x-1)/2)*imgData.pxSzT (c_m-1+(n_px_x-1)/2)*imgData.pxSzT], ...
            get(ha_img,'YLim'), 'Parent',ha_img, ...
            'LineWidth',1, 'Color','k', 'LineStyle','-');
        set(ha_img, 'FontSize',fontSzNum)
        old_img_cLims = get(ha_img, 'Clim');
        % in case taht in image of events there is also a part
        % of another event with much higher amplitude
        set(ha_img, 'Clim', [old_img_cLims(1) ...
            prctile(evnt.axData{1}.imgE(evnt.axData{1}.imgE_m),99,'all')])
        title(ha_img, evnt.axDesc{1}.fig_title_txt, 'FontSize',fontSzT)
        xlabel(ha_img,'t (ms)', 'FontSize',fontSzL)
        ylabel(ha_img,'x (\mum)', 'FontSize',fontSzL)
end

% x-profile axes
ha_xProf.YAxis.Direction = "reverse";
line(evnt.axData{1}.x_event_prof, evnt.axData{1}.x, ...
    'Parent',ha_xProf, ...
    'LineWidth',1, 'Color','k', 'Display','data')
line(evnt.axData{1}.x_prof_fit, evnt.axData{1}.x_ups, ...
    'Parent',ha_xProf, ...
    'LineWidth',2, 'Color','r', 'Display',evnt.axData{1}.x_prof_fit_txt)
set(ha_xProf,'FontSize',fontSzNum)
title(ha_xProf, evnt.axDesc{1}.ax_profX_title_txt, 'FontSize',fontSzT)
xlabel(ha_xProf, hObjs_mainFig.ax_prof.YLabel.String, 'FontSize',fontSzL)
ylabel(ha_xProf, 'x (\mum)', 'FontSize',fontSzL)
ylim(ha_xProf, [min(evnt.axData{1}.x_ups) max(evnt.axData{1}.x_ups)])
xlim(ha_xProf, getAxisLimits( ...
    [evnt.axData{1}.x_event_prof(:);evnt.axData{1}.x_prof_fit(:)], 5))
line([evnt.axData{1}.half_max_x evnt.axData{1}.half_max_x], ...
    [evnt.axData{1}.half_max_x_1 evnt.axData{1}.half_max_x_2], ...
    'Parent',ha_xProf,...
    'LineWidth',2, 'Color','b', 'Display','FWHM')
hl_2 = legend(ha_xProf,'show');
hl_2.Location = 'best';
hl_2.FontSize = fontSzLegend;

% t-profile axes
line(evnt.axData{1}.t, evnt.axData{1}.t_event_prof, 'Parent',ha_tProf, ...
    'LineWidth',1, 'Color','k', 'Display','data')
line(evnt.axData{1}.t_ups, evnt.axData{1}.t_prof_fit, 'Parent',ha_tProf, ...
    'LineWidth',2, 'Color','r', 'Display',evnt.axData{1}.t_prof_fit_txt)
% show exponential fits
if ~isempty(evnt.axData{1}.t_rise)
    line(evnt.axData{1}.t_rise, evnt.axData{1}.t_profR_fit, ...
        'Parent',ha_tProf, ...
        'LineWidth',2, 'Color','c', ...
        'LineStyle',':', 'Display','exp rise fit')
    line(evnt.axData{1}.t_decay, evnt.axData{1}.t_profD_fit, ...
        'Parent',ha_tProf, ...
        'LineWidth',2, 'Color','m', ...
        'LineStyle',':', 'Display','exp decay fit')
end
set(ha_tProf, 'FontSize',fontSzNum)
title(ha_tProf, evnt.axDesc{1}.ax_profT_title_txt, 'FontSize',fontSzT)
xlabel(ha_tProf, 't (ms)', 'FontSize',fontSzL)
ylabel(ha_tProf, hObjs_mainFig.ax_prof.YLabel.String, 'FontSize',fontSzL)
xlim(ha_tProf, [min(evnt.axData{1}.t_ups) max(evnt.axData{1}.t_ups)])
ylim(ha_tProf, getAxisLimits( ...
    [evnt.axData{1}.t_event_prof(:);evnt.axData{1}.t_prof_fit(:)], 5))
line([evnt.axData{1}.half_max_t_1 evnt.axData{1}.half_max_t_2], ...
    [evnt.axData{1}.half_max_t evnt.axData{1}.half_max_t], 'Parent',ha_tProf,...
    'LineWidth',2, 'Color','b', 'Display','FDHM')
line([evnt.axData{1}.t_max evnt.axData{1}.t_max], ...
    [evnt.axData{1}.bs_t evnt.axData{1}.v_max], 'Parent',ha_tProf, ...
    'LineWidth',2, 'Color','g', 'Display','amplitude')
line([evnt.axData{1}.t0 evnt.axData{1}.t_max], ...
    [evnt.axData{1}.bs_t evnt.axData{1}.bs_t], 'Parent',ha_tProf, ...
    'LineWidth',2, 'Color','m', 'Display','TTP')
hl_3 = legend(ha_tProf, 'show');
hl_3.Location = 'best';
hl_3.FontSize = fontSzLegend;

% show parameters of analyzed event
h_txt_params =  hObjs_evntsBrowser.h_txt_params;
h_txt_params.String = {
    sprintf('calcMethod: %s', evnt.evntParams{1}.calcMethod{1}), ...
    [sprintf('amplitude = %0.2f',evnt.evntParams{1}.amplitude), ...
    ' (',char(916),sprintf('F/F0)')],...
    sprintf('TTP = %0.2f (ms)',evnt.evntParams{1}.TTP), ...
    sprintf('FDHM = %0.2f (ms)',evnt.evntParams{1}.FDHM), ...
    [sprintf('FWHM = %0.2f ',evnt.evntParams{1}.FWHM), ...
    '(',char(181),'m)'],...
    [sprintf('sparkMass = %0.2f ',evnt.evntParams{1}.sparkMass), ...
    '(',char(916),'F/F0*',char(181),'m^3)'], ...
    sprintf('accepted event = %d ',evnt.maskOfAcceptedSparks)};

h_txt_params.HorizontalAlignment = 'left';
h_txt_params.FontSize = 0.075;
h_txt_params.FontWeight = 'normal';

end