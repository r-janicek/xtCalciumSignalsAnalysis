function bsMask_slider(h_sld,~,fitFig)


% get slider value
sldVal = round(h_sld.Value);
if sldVal==100, sldVal=99; end
if sldVal==0, sldVal=1; end

hObjsFit = getappdata(fitFig,'hObjsFit');

% get profile data
selectedProf = getappdata(fitFig,'selectedProf');

% calculate break points for detrending
N_intervals = round(selectedProf.t(end)/1000);
if N_intervals<2, N_intervals=2; end
bp = (1:1:N_intervals-1).*(round(numel(selectedProf.t)/N_intervals));
bp = [bp,numel(selectedProf.t)];
y_d = detrend(selectedProf.y,'linear',bp);

% calculate and save new mask
% mBs = selectedProf.y <= prctile( selectedProf.y, sldVal );
% calculate baseline mask from detrended data, better estimate
mBs = y_d <= prctile( y_d, sldVal );

if strcmp(hObjsFit.popUpMenuBs.String{hObjsFit.popUpMenuBs.Value},'spline')
   
    mBs(1) = true;
    mBs(end) = true;
    
end
        
selectedProf.baselineM = mBs(:);
setappdata(fitFig,'selectedProf',selectedProf);


% do fitting of baseline
fitBaseline([],[],fitFig)

end

