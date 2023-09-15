function selectBaselineFcn(hO,~,mainFig)

% get data
imgData = getappdata(mainFig,'imgData');
hObjs = getappdata(mainFig,'hObjs');

profLineObj = findall(hObjs.ax_prof,'Tag','wholeCellProfile');
prof = profLineObj.YData;
t =  profLineObj.XData;

type = hO.String{hO.Value};

% do fitting of baseline
bsFit = fitBaseline(t,prof,type,imgData.baselineM,1,hObjs.ax_prof);

% save new fit
imgData.bsFit = bsFit;
setappdata(mainFig,'imgData',imgData);

end

