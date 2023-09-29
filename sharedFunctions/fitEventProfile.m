function outFitEvent = fitEventProfile(t,y,analysisFig,axFit,axRes,tauD_regionWise,waveSpeed)

% ax = axes to plot 
hObjsA = getappdata(analysisFig,'hObjsA');
N_events = str2double(hObjsA.h_edit_Npeaks.String);
pxSzT = mean(diff(t));
t_ups = (min(t):pxSzT/10:max(t));
keyboard
paramFit1 = str2double(hObjsA.h_edit_paramFit1.String);
paramFit2 = str2double(hObjsA.h_edit_paramFit2.String);
popUpFitFun = hObjsA.popUpMenuEventsFcn; 
model = popUpFitFun.String{popUpFitFun.Value};
popUpBsFun = hObjsA.popUpMenuBaselineFcn; 
modelBs = popUpBsFun.String{popUpBsFun.Value};
pieceWise = hObjsA.check_pieceWise.Value; 
if pieceWise
    pieceWise='yes'; 
else
    pieceWise='no';
end

if hObjsA.check_signOfPeak.Value
    signOfPeak = '+'; 
else
    signOfPeak = '-'; 
end

% smooth profile
sp = (str2double(hObjsA.h_smooth_edit.String)/pxSzT)/numel(y);  % span <0,1>
ys = smooth(y,sp,'loess');

% get initial guess for fitting
[pks,locs,w,p] = findpeaks(ys,t,'NPeaks',N_events,'SortStr','descend');

%[vMax,posMax] = max(smooth(y,5));
t0 = locs(:) - w(:)./2;
tpeak = locs(:);
F0 = zeros(size(t0));
tauR = w(:)./3;
tauD = w(:);
A = pks(:);

switch model
    
    case 'spline'
        % coefficients: [t0, F01, spline coefficients]
        p0 = [t0,F0];
        
    case 'expModGauss'
        % coefficients: [A,m,sd,tau,F0]
        p0 = [A,tpeak,w(:)./2,tauD,F0];
        if strcmp(pieceWise,'yes')
            errordlg('piecewise is not supported for expModGauss')
            return
        end
        
    case 'doubleBoltzmann'
        % coefficients: [F0, A, xH1, k1, xH2, k2]
        p0 = [F0,A,tpeak-w(:)./2,w(:)./2,tpeak+w(:)./2,-w(:)./2];
        if strcmp(pieceWise,'yes')
            errordlg('piecewise is not supported for doubleBoltzmann')
            return
        end
        
    case 'Asym2Sigmoid'
        % coefficients:[F0, A, xc, w1, w2, w3]
        % Meanings: y0 = offset, xc = center, A = amplitude, w1 = full width of half maximum,
        % w2 = variance of low-energy side, w3 = variance of high-energy side.
        % Lower Bounds: w1 > 0.0, w2 > 0.0, w3 > 0.0
        p0 = [F0,A,tpeak,w(:),w(:)/2,w(:)/2];
        if strcmp(pieceWise,'yes')
            errordlg('piecewise is not supported for Asym2Sigmoid')
            return
        end   
        
    case '1expR1expD'
        % coefficients: [t0, F01, tauR, A, t1, tauD, F02]
        p0 = [t0,F0,tauR,A,tpeak,tauD,F0];
        
    case 'CaSpikeFun'
        %coefficients: [t0, F0, FM, tA, tI, FI]
        p0 = [t0,F0,A,tauR,tauD,ones(size(t0))];
        
end
%keyboard

% fit profile
% delete previous fits and lines
delete(findobj(axFit,'Tag','fitProfile'))
delete(findobj(axRes,'Tag','fitResiduals'))
delete(findobj(axFit,'Tag','params'))

if N_events>1
    % multiple peaks fitting
    
    fittedEvent = multipleEventsFitFun(t,y,ys,N_events,[],[],modelBs,...
        [],model,t_ups,3,signOfPeak,[],[],axFit,0);

else
    % single peak fitting
    fittedEvent = fitOneWholeEvent(p0,t,t_ups,...
        y,pieceWise,model,paramFit1,paramFit2);

end

% get parameters of fitted peaks in profile
calcParams = calcParametersOfEventsFromWholeProfileFit(t,fittedEvent);
outFitEvent = fittedEvent;
outFitEvent.calcParamsFromFit = calcParams;
 
%keyboard

% plot result
delete(findall(axFit,'Tag','fitProfile'))
try
line(t,outFitEvent.yNorm,...
    'Parent',axFit,'Color','k','LineStyle','-','LineWidth',1,...
    'Tag','fitProfile');
line(t,outFitEvent.wholeFit,...
    'Parent',axFit,'Color','r','LineStyle','-','LineWidth',2,...
    'Tag','fitProfile');
line(t,outFitEvent.ysNorm,...
    'Parent',axFit,'Color','b','LineStyle','-','LineWidth',1,...
    'Tag','fitProfile');
line([calcParams{2:end,3}],[calcParams{2:end,7}],...
    'Parent',axFit,'Color','c','LineStyle','none','LineWidth',1,...
    'Marker','.','MarkerSize',30,...
    'Tag','fitProfile');

plot(t,outFitEvent.ysNorm(:)-outFitEvent.wholeFit(:),...
    'Parent',axRes,'Color','k','LineStyle','-','LineWidth',1,...
    'Tag','fitResiduals');

% set axes limits
axFit.XLim = [min(t) max(t)];
axFit.YLim = [min(outFitEvent.yNorm)*0.95 max(outFitEvent.yNorm)*1.05];

axRes.XLim = [min(t) max(t)];
%axRes.YLim = [ ];
catch
end
bsFit = fittedEvent.t_ups.baselineFit;
% show phenomenological parameters
for i=2:size(calcParams,1)
    try
        bs = bsFit(round(calcParams{i,1}/pxSzT));
    catch
        bs = bsFit(1);
    end
    line([calcParams{i,1} calcParams{i,3}],[calcParams{i,2}+bs calcParams{i,2}+bs],...
        'Parent',axFit,'Color','b','LineStyle','-','LineWidth',2,...
        'Tag','params')
    line([calcParams{i,4} calcParams{i,5}],[calcParams{i,6}+bs calcParams{i,6}+bs],...
        'Parent',axFit,'Color','b','LineStyle','-','LineWidth',2,...
        'Tag','params')
    line([calcParams{i,3} calcParams{i,3}],[calcParams{i,2}+bs calcParams{i,7}+bs],...
        'Parent',axFit,'Color','b','LineStyle','-','LineWidth',2,...
        'Tag','params')
    line(calcParams{i,14}(:,1), calcParams{i,14}(:,2),...
        'Parent',axFit,'Color','b','LineStyle',':','LineWidth',2,...
        'Tag','params')   
end

if size(calcParams,1)==2
    text(axFit,t(5),max(y),{[sprintf('Amplitude: %0.2f',calcParams{2,7}),'  \DeltaF/F_0'];...
        sprintf('TTP: %0.2f ms',calcParams{2,8});...
        sprintf('FDHM: %0.2f ms',calcParams{2,9});...
        sprintf('tauD fit: %0.2f ms',calcParams{2,11});...
        sprintf('tauD median pxWise: %0.2f ms',median(tauD_regionWise));...
        [sprintf('wave speed: %0.2f',waveSpeed),' \mum/s']},...
        'VerticalAlignment','top','Tag','params',...
        'FontSize',14)
end


end










