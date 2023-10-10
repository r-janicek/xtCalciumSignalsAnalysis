function fitBaseline(~, ~, fitFig)

set(fitFig,'Pointer','watch')
drawnow

hObjsFit = getappdata(fitFig,'hObjsFit');

% get selected profile data
prof = getappdata(fitFig,'selectedProf');

switch hObjsFit.popUpMenuBs.String{hObjsFit.popUpMenuBs.Value}
    
    case 'spline'
        n_knots = str2double(hObjsFit.h_edit_paramFitBs1.String);
        splOrd = str2double(hObjsFit.h_edit_paramFitBs2.String);  % spline order
        
        spl = spap2(n_knots,splOrd,prof.t(:),prof.y(:),double(prof.baselineM(:)));
        % newknt for a possibly better knot distribution
        knots = newknt(spl);
        % least-squares approximation
        spl = spap2(knots, splOrd, prof.t, prof.y, double(prof.baselineM(:)));
        bsFit = fnval(spl,prof.t);
        
    
    case 'SmoothingSpline'
        s_span = str2double(hObjsFit.h_edit_paramFitBs1.String);
        ft = fittype('smoothingspline');
        opts = fitoptions('Method','SmoothingSpline',...
            'SmoothingParam',s_span,...
            'Normalize','on',...
            'Exclude',~prof.baselineM);
        % fit
        [f_bs, ~,~] = fit(prof.t(:), prof.y(:), ft, opts);
        bsFit = feval(f_bs,prof.t);
     
        
    case 'polynomial' 
        order = hObjsFit.h_edit_paramFitBs1.String;
        ft = fittype(['poly',order]);
        opts = fitoptions('Method','LinearLeastSquares',...
            'Normalize','on',...
            'Robust','Bisquare',...
            'Exclude',~prof.baselineM);
        
        % fit
        [f_bs, ~,~] = fit(prof.t(:), prof.y(:), ft, opts);
        bsFit = feval(f_bs,prof.t);
        
            
    case 'stretchedExp'
        fitFun = @(p,x) p(1).*exp( -(x./p(2)).^p(3) ) + p(4);
        fitFunSum = @(p,t,y) sum( (y-fitFun(p,t)).^2 );
        % do fit
        opt = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
            'MaxIter',1000,'MaxFunEvals',3000);
        coef = fmincon(@(p)fitFunSum(p,prof.t(prof.baselineM),prof.y(prof.baselineM)),...
            [100 prof.t(end)/10 0.5 median(prof.y(prof.baselineM))],[],[],[],[],[-inf -inf  0 -inf],[inf inf 1 inf],[],opt);
%         mean(prof.y(find(prof.baselineM,100,'last')))
%            figure
%            plot(prof.t,prof.y)
%            hold on
%            plot(prof.t,fitFun(coef,prof.t),'r')
%            plot(prof.t,fitFun([1 100 0.25 median(prof.y(prof.baselineM))],prof.t),'g')
                
        % baseline fit
        bsFit = fitFun(coef,prof.t);
        
        
    otherwise
        ft = fittype(hObjsFit.popUpMenuBs.String{hObjsFit.popUpMenuBs.Value});
        opts = fitoptions('Method','NonlinearLeastSquares',...
            'Normalize','on',...
            'Robust','Bisquare',...
            'Exclude',~prof.baselineM);
        
        % fit
        [f_bs, ~,~] = fit(prof.t(:), prof.y(:), ft, opts);
        bsFit = feval(f_bs,prof.t);
        
end

% delete previous fit
delete(findobj(hObjsFit.ax_fit,'Tag','baselineFit'))
% new fit
line(prof.t,bsFit,'Parent',hObjsFit.ax_fit,'Color','b','LineStyle','-','LineWidth',2,...
    'Tag','baselineFit');

% delete previous mask
hl_m_Old = findobj(hObjsFit.ax_fit,'Type','Line','-regexp','Tag','Mask');
delete(hl_m_Old)
% new mask
hl_m = line(prof.t(prof.baselineM),prof.y(prof.baselineM),'Parent',hObjsFit.ax_fit,'Color','g',...
                     'LineStyle','none','Marker','.','MarkerSize',20,'LineWidth',1,'Tag','baselineMask');
uistack(hl_m, 'bottom')

% save baseline fit
prof.baselineFit = bsFit;
setappdata(fitFig,'selectedProf',prof);

% update baseline slider
bsSldVal = round(sum(prof.baselineM)/numel(prof.baselineM)*100);
if bsSldVal==100, bsSldVal=99; end
if bsSldVal==0, bsSldVal=1; end
hObjsFit.h_txt_sld_Bs.String = sprintf('mask of baseline as %d percentile',bsSldVal);
hObjsFit.sld_Bs.Value = bsSldVal;


set(fitFig,'Pointer','arrow')
drawnow

end

