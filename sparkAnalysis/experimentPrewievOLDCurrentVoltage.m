function experimentPrewievOLDCurrentVoltage(~,~,mainFig)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% pridat save current, voltage,... do struktury,
% save parametrov do tabulky
keyboard
% get data 
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');

if isempty(imgData)
    openImg([],[],mainFig)
end

if isfield(imgData,'filePath')
    
    [FileNames, FilePath] = uigetfile({'*.bwav'},'Load voltage & current related to image',...
        imgData.filePath,'Multiselect','on');
else
    [FileNames, FilePath] = uigetfile({'*.bwav'},'Load voltage & current related to image',...
        'c:\Documents and Settings\Rado Janicek\Desktop\',...
        'Multiselect','on');
end

FilePath_U = [FilePath,FileNames{1,2}];
FilePath_C = [FilePath,FileNames{1,1}];

voltage = IBWread(FilePath_U);
curr = IBWread(FilePath_C);

img = imgData.wholeImgFluoXT;
px_t = imgData.pxSzT;
xt_cur = (0:1:numel(curr.y)-1).*curr.dx;

scRes = get(0,'ScreenSize');

% plot output
whImg_fig = figure;
set(whImg_fig,'Position',[1 1 scRes(3) scRes(4)],'Name','Experiment prewiev','Tag','Experiment_prewiev')

vol_ax = axes('Parent',whImg_fig,'Units','normalized','Position',[0.03 0.89 0.94 0.10]);
curr_ax = axes('Parent',whImg_fig,'Units','normalized','Position',[0.03 0.73 0.94 0.15]);
whImg_ax = axes('Parent',whImg_fig,'Units','normalized','Position',[0.03 0.47 0.94 0.25]);
prof_ax = axes('Parent',whImg_fig,'Units','normalized','Position',[0.03 0.26 0.94 0.2]);
text_ax = axes('Parent',whImg_fig,'Units','normalized','Position',[0.22 0.01 0.1 0.15]);
set(text_ax,'Visible','off')

linkprop([vol_ax,curr_ax,whImg_ax,prof_ax],'XLim');

electroPhys = struct('xt_cur',xt_cur,...
                     'voltage',voltage.y,...
                     'current',curr.y);  

setappdata(mainFig,'electroPhys',electroPhys)
                
plotExperimentPrewiev([],[],mainFig);

h_pb_delay = uicontrol('Style', 'pushbutton',...
    'String','set delay','FontUnits','normalized','FontSize',0.2,'FontWeight','bold',...
    'Parent',whImg_fig,'Units','normalized','Position', [0.8 0.11 0.1 0.1],...
    'Callback',{@plotExperimentPrewiev,mainFig},'Enable','on');

% table of selected profiles
h_table = uitable('Data',[{'Cm'; 'Rm'; 'Ra'; 'tau'}, cellstr( num2str(nan(4,1)) ), {'pF'; 'MOhms'; 'MOhms'; 'us'}],...
    'ColumnEditable',[false false false],...
    'ColumnName',{'<HTML> <font size="5"> <b> parameter </b> </HTML>', ...
    '<HTML> <font size="5"> <b> ---value--- </b> </HTML>', '<HTML> <font size="5"> <b> ---units--- </b> </HTML>'},...
    'ColumnWidth','auto','Parent',whImg_fig,...
    'Units','normalized','Position',[0.225 0.05 0.175 0.16],...
    'FontUnits','normalized','FontSize',0.10,...
    'CellEditCallback','');

h_pb_cellParam = uicontrol('Style', 'pushbutton',...
    'String','ROI for cell params','FontUnits','normalized','FontSize',0.2,'FontWeight','bold',...
    'Parent',whImg_fig,'Units','normalized','Position', [0.1 0.11 0.1 0.1],...
    'Callback',{@calcCellParams,mainFig,h_table},'Enable','on');

check_pulseDirection = uicontrol('Style','checkbox','Parent',whImg_fig,...
    'FontUnits','normalized','Value',1,'Units','normalized','FontSize',0.3,...
    'String','negative pulses','Position', [0.1 0.05 0.4 0.05]);

h_pb_close = uicontrol('Style', 'pushbutton',...
    'String','close & save data','FontUnits','normalized','FontSize',0.2,'FontWeight','bold',...
    'Parent',whImg_fig,'Units','normalized','Position', [0.425 0.11 0.1 0.1],...
    'Callback',{@closeFun,mainFig,whImg_fig},'Enable','on');



%%%%%%%%%%%%
    function plotExperimentPrewiev(~,~,mainFig)
        
        %         keyboard
        %         pridat automaticke hladanie delay
        %         startOfFirstPulse = find(round(voltage.y./5).*5 == max(round(voltage.y./5).*5),1,'first');
        %         findpeaks( mean(img,1)   )
        %
        %         figure
        %         plot(round(voltage.y./5).*5)
        
        % get a delay
        answer = str2double(inputdlg1({'Insert delay between image and current (ms):','caffeine (0-no, 1-yes)'},...
            'delay',1,{'9500','0'}));
        
        delayImgToCurr = answer(1,1);
        caffeine = logical(answer(2,1));
        
        xt_img = linspace(0,(size(img,2)-1)*px_t,size(img,2)) + delayImgToCurr;
        %plot voltage
        plot(xt_cur,voltage.y,'Parent',vol_ax,'Color','k');
        ylabel(vol_ax,'voltage (mV)')
        set(vol_ax,'XTick',[],'FontSize',16)
        set(vol_ax,'XLim',[0 max(xt_img(end),xt_cur(end))],'YLim',[min(voltage.y)-5 max(voltage.y)+5])
        
        if delayImgToCurr/curr.dx<1
            curr_img = curr.y(1:end);
            volt_img = voltage.y(1:end);
            
        else
            curr_img = curr.y(delayImgToCurr/curr.dx:end);
            volt_img = voltage.y(delayImgToCurr/curr.dx:end);
        end
        
        if caffeine
            
            nPulse = 0;
            ICa_ampl = 0;
            
        else
            
            % find mask of voltage waveform
            volt_img_mask = volt_img;
            volt_img_mask = round(volt_img_mask./5).*5;
            
            % find possible start of test pulses (>=-35mV)
            startOfTrainOfPulses = find(volt_img_mask>=-35,1,'first');
            
            % add baseline before first pulse
            startOfTrainOfPulses = find(diff(volt_img_mask(1:startOfTrainOfPulses-1))~=0,1,'last')+1;
            
            if isempty(startOfTrainOfPulses), startOfTrainOfPulses=1; end
            
            % find end of pulses
            endOfTrainOfPulses = find(volt_img_mask>=-35,1,'last')+1;
            
            % find start of pulses
            startOfIndividualPulses = find(diff(volt_img_mask)>0);
            startOfIndividualPulses(startOfIndividualPulses < startOfTrainOfPulses)=[];
            startOfIndividualPulses(startOfIndividualPulses > endOfTrainOfPulses+1)=[];
            
            % find end of pulses
            endOfIndividualPulses = find(diff(volt_img_mask)<0);
            endOfIndividualPulses(endOfIndividualPulses < startOfTrainOfPulses)=[];
            endOfIndividualPulses(endOfIndividualPulses > endOfTrainOfPulses+1)=[];
            
            % volt_img_mask = volt_img;
            % volt_img_mask(volt_img_mask < -40) = -40;
            % [idx,C] = kmeans(volt_img_mask(volt_img_mask>-45),2);
            
        end
        
        %         [~,p] = max(C);
        %
        %         if p==1
        %
        %             diff_volt_img_pos = find(diff(idx)<0);
        %             diff_volt_img_neg = find(diff(idx)>0);
        %
        %         else
        %             diff_volt_img_pos = find(diff(idx)>0);
        %             diff_volt_img_neg = find(diff(idx)<0);
        %         end
        %         keyboard
        
        %         figure
        %         plot(volt_img_mask)
        %
        % calculate ICa aplitude, position and baseline for individual
        % pulses
        if caffeine
            
            YCurrLlimits = [min(curr_img),max(curr_img)];
            
        else
            
            for i=1:numel(startOfIndividualPulses)
                
                if i==1
                    
                    curr_pulse = curr_img(startOfIndividualPulses(i):endOfIndividualPulses(i));
                    curr_baseline = curr_img(startOfTrainOfPulses:startOfIndividualPulses(i));                                                                               
                else
                    
                    curr_pulse = curr_img(startOfIndividualPulses(i):endOfIndividualPulses(i));
                    curr_baseline = curr_img( round((startOfIndividualPulses(i) - endOfIndividualPulses(i-1))/2) + endOfIndividualPulses(i-1) : startOfIndividualPulses(i));
                end
                
                voltageStep(i,1) = max(volt_img_mask(startOfIndividualPulses(i):endOfIndividualPulses(i)));
                                
                ICa_baseline(i,1) = median(curr_baseline);
                
                [ICa_peak(i,1),ICa_peak_pos(i,1)] = min(smooth(curr_pulse,ceil(1/curr.dx)));
                
                ICa_peak_pos(i,1) = ICa_peak_pos(i,1)+startOfIndividualPulses(i);
                
                nPulse(i,1) = {sprintf('%d pulse',i)};
                ICa_ampl(i,1) = sqrt((ICa_peak(i,1) - ICa_baseline(i,1))^2);
                
                YCurrLlimits(i,:) = [min(min(curr_baseline),min(curr_pulse)),...
                    max(max(curr_baseline),max(curr_pulse(round(end/2):end-2)))];
                
                clearvars curr_pulse curr_baseline
            end
        end
        %plot voltage, current and image
        plot(xt_cur,curr.y,'Parent',curr_ax,'Color','r');
        ylabel(curr_ax,'current (pA)')
        set(curr_ax,'XTick',[],'FontSize',16)
        
        minIcaYlim = min(YCurrLlimits(:,1))-100;
        maxIcaYlim = max(YCurrLlimits(:,2))+100;
        
        set(curr_ax,'XLim',[0 max(xt_img(end),xt_cur(end))],...
            'YLim',[minIcaYlim, maxIcaYlim])
        
        if ~caffeine
            for i=1:numel(ICa_peak)
                
                pos_x = ICa_peak_pos(i,1)*curr.dx + delayImgToCurr-50;
                
                line([pos_x pos_x],[ICa_baseline(i,1) ICa_peak(i,1)],'Parent',curr_ax,'LineWidth',2,'Color','k',...
                    'LineStyle',':');
                
            end
        end
        image(img,'CDataMapping','scaled','Parent',whImg_ax,'XData',[xt_img(1) xt_img(end)]);
        set(whImg_ax,'XLim',[0 max(xt_img(end),xt_cur(end))])
        colormap(parula(256))
        set(whImg_ax,'FontSize',16,'XTick',[])
        ylabel(whImg_ax,'x (pixels)')
        
        plot(xt_img,mean(img,1),'Parent',prof_ax,'Color','b');
        ylabel('current (pA)')
        set(prof_ax,'FontSize',16)
        set(prof_ax,'XLim',[0 max(xt_img(end),xt_cur(end))])
        ylabel(prof_ax,'mean fluorescence')
        xlabel(prof_ax,'t (ms)')
        
        %plot position of TPP
        try
            crop_s_t = imgData.crop_s_t;
            if isempty(crop_s_t), crop_s_t=0; end
            
            s_TPP = imgData.s_TPP;
            e_TPP = imgData.e_TPP;
            
            for i=1:numel(s_TPP)
                
                line([xt_img(s_TPP(i)) xt_img(s_TPP(i))]+crop_s_t,get(whImg_ax,'YLim'),'Parent',whImg_ax,'LineWidth',2,'Color','r',...
                    'LineStyle',':');
                line([xt_img(e_TPP(i)) xt_img(e_TPP(i))]+crop_s_t,get(whImg_ax,'YLim'),'Parent',whImg_ax,'LineWidth',2,'Color','r',...
                    'LineStyle',':');
            end
        catch
        end
        
        % save Ica values
        electroPhys = getappdata(mainFig,'electroPhys');
        electroPhys.currentAmpl = [[{'voltage pulse #'};nPulse],[{'V (mV)'};num2cell(voltageStep)],[{'Ica (pA)'};num2cell(ICa_ampl)]];
        electroPhys.delayImgToCurr = delayImgToCurr;
        electroPhys.IcaLimits = [minIcaYlim, maxIcaYlim];
        
        setappdata(mainFig,'electroPhys',electroPhys)
        
        linkprop([vol_ax,curr_ax,whImg_ax,prof_ax],'XLim');
        
    end


%%%%%%%%%%%%%%
    function calcCellParams(h_O,~,mainFig,h_table)
        
        % fit: y = Iss + B*exp(-x/tau)
        % Iss = U/(Rm + Rs)
        % tau = Cm*Rs*Rm/(Rs + Rm)
        % Rs = U/(Ipeak-mean(I1))
        
        % delay of current due to filter (lowpass Bessel,10 kHz)
        % measured on model cell in ms, for ~ 10 mV pulse
        % and model cell parameters Rm = 500 mOhms, Rs = 10 mOhms, Cp = 33pF
        %Rs = 10;
        %f_delay(i,1) = -coeff(3)*log((1000*U_step + abs(coeff(1))*Rs + mean(I1)*Rs)/(coeff(2)*Rs));
        
        peakOrientation = (-1)^check_pulseDirection.Value;
        st = h_O.String;
        
        switch st
            
            case 'ROI for cell params'
                
                xL = get(curr_ax,'XLim');
                yL = get(curr_ax,'YLim');
                
                h_rect_params = imrect(curr_ax,[xL(1) yL(1) xL(2)/10 abs(yL(1))+abs(yL(2))]);
                setColor(h_rect_params,'b')
                % constrains for ROI
                fcn = makeConstrainToRectFcn('imrect',xL,yL);
                setPositionConstraintFcn(h_rect_params,fcn);
                setResizable(h_rect_params,1);
                
                setappdata(whImg_fig,'h_rect_params',h_rect_params)
                
                h_O.String = 'get cell params';
                
            case 'get cell params'
                
                h_rect_params = getappdata(whImg_fig,'h_rect_params');
                try
                    pos = getPosition(h_rect_params);
                    pos = round([pos(1)/curr.dx pos(2) pos(3)/curr.dx pos(4)]); %convert to points
                    
                    pos(pos<1) = 1;
                    
                    V = voltage.y;
                    I = curr.y;
                    
                    % during pulse
                    V_p = V(pos(1):pos(1)+pos(3));
                    I_p = I(pos(1):pos(1)+pos(3));
                    
                    %shift I_p to positive values
                    switch peakOrientation
                        case 1
                            I_p = I_p + abs(min(I_p));
                        case -1
                            I_p = abs(I_p);
                    end
                    
                    % find two groups baseline and peak from voltage waveform
                    [idx,C] = kmeans(V_p,2);
                    
                    % find positions of voltage changing
                    idx_change = find(diff(idx)~=0);
                    
                    % voltage step
                    U_step = abs(abs(C(1)) - abs(C(2)));
                    %steady state current before first capacitance peak
                    I1 = I_p(1:idx_change(1));
                    
                    for i = 1:numel(idx_change)/2
                        
                        I_pp = I_p(idx_change(2*i-1):idx_change(2*i));
                        
                        % remove part before peak
                        [v_p,p_p] = max(I_pp);
                        I_pp = I_pp(p_p:end);
                        t = (0:1:numel(I_pp)-1).*curr.dx;
                        
                        % start point
                        [~,tau_p] = min(abs(v_p*0.5 - I_pp));
                        x0 = [I_pp(end) v_p t(tau_p)];
                        fo = fitoptions('Method','NonlinearLeastSquares','Robust','off','StartPoint',x0);
                        ft = fittype('Iss + B*exp(-t/tau)','Coefficient',{'Iss', 'B', 'tau'},'Independent','t');
                        
                        f_ = fit(t(:),I_pp(:),ft,fo);
                        I_fit = feval(f_,t);
                        coeff = coeffvalues(f_);
                        Iss = coeff(1);
                        
                        deltaIss = abs(Iss - mean(I1));
                        Ipeak = abs(max(I_fit) - mean(I1));
                        
                        %charge in pA*ms
                        %keyboard % check for second term in equation, it should be correction factor for calculation of capacitance
                        Q = sum((I_fit-Iss).*curr.dx) + Iss*coeff(3);
                        
                        % %capacitance of the cell in pF
                        Cm(i,1) = Q/U_step;
                        
                        %access resistance in MOhms
                        Ra(i,1) = (U_step/Ipeak)*1000;
                        
                        %cell resistance in MOhms
                        Rm(i,1) = (U_step/deltaIss)*1000 - Ra(i,1);
                        
                        %time constant in us
                        tau(i,1) = coeff(3)*1000;
                        
                        % % more accurate ???
                        %Cm(i,1) = (tau(i,1)*(Ra(i,1)+Rm(i,1)))/(Ra(i,1)*Rm(i,1));
                        
                        ff(i) = figure;
                        plot(t,I_pp,'ko')
                        hold on
                        plot(t,I_fit,'r')
                        
                        clearvars I_pp t I_fit
                        
                    end
                    
                    Cm(  Cm<=0  | isnan(Cm)  )=[];
                    Rm(  Rm<=0  | isnan(Rm)  )=[];
                    Ra(  Ra<=0  | isnan(Ra)  )=[];
                    tau( tau<=0 | isnan(tau) )=[];
                    
                catch
                    Cm = nan;
                    Rm = nan;
                    Ra = nan;
                    tau = nan;
                end
                
                cellParameters = table(mean(Cm),mean(Rm),mean(Ra),mean(tau),...
                    'VariableNames',{'Cm', 'Rm', 'Ra', 'tau'});
                cellParameters.Properties.VariableUnits = {'pF' 'MOhms' 'MOhms' 'us'};
                
                electroPhys = getappdata(mainFig,'electroPhys');
                electroPhys.cellParameters = cellParameters;
                setappdata(mainFig,'electroPhys',electroPhys)
                
                h_table.Data = [{'Cm'; 'Rm'; 'Ra'; 'tau'},...
                    cellstr( num2str([mean(Cm);mean(Rm);mean(Ra);mean(tau)]) ),...
                    {'pF'; 'MOhms'; 'MOhms'; 'us'}];
                
                if exist('ff'), close(ff); end
                delete(h_rect_params)
                h_O.String = 'ROI for cell params';
                
        end
        
    end


%%%%%%%%%%%%%%
    function closeFun(~,~,mainFig,whImg_fig)
        
        electroPhys = getappdata(mainFig,'electroPhys');
        
        if ~isfield(electroPhys,'cellParameters')
            
            cellParameters = table(nan,nan,nan,nan,...
                'VariableNames',{'Cm', 'Rm', 'Ra', 'tau'});
            cellParameters.Properties.VariableUnits = {'pF' 'MOhms' 'MOhms' 'us'};
            
            electroPhys.cellParameters = cellParameters;
            setappdata(mainFig,'electroPhys',electroPhys)
        end
        
        close(whImg_fig)
    end


end

