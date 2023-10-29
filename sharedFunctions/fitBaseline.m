function bsFit = fitBaseline(t, y, type, mBs, plotResult, ax, params)

% t = time
% y = profile
% type = type of function
% mBs = mask of baseline
% plotResult = 1 or 0, plot results in axes ax
% params = parameters of fit function 
if strcmp(type, 'polynomial')
    type = ['poly', num2str(params(1))];
end

switch type
       
    case {'poly1','poly2','poly3','poly4',...
            'poly5','poly6','poly7','poly8','poly9'}

        ft = fittype(type);
        opts = fitoptions('Method','LinearLeastSquares',...
            'Normalize','on',...
            'Robust','Bisquare',...
            'Exclude',~mBs);
        % fit
        [f_bs, ~,~] = fit(t(:), y(:), ft, opts);
        bsFit = feval(f_bs, t);
        
            
    case 'stretchedExp'
        fitFun = @(p,x) p(1).*exp( -(x./p(2)).^p(3) ) + p(4);
        fitFunSum = @(p,t,y) sum( (y-fitFun(p,t)).^2 );
        % do fit
        opt = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
            'MaxIter',1000,'MaxFunEvals',3000);
        coef = fmincon(@(p)fitFunSum(p, t(mBs), y(mBs)),...
            [100 t(end)/10 0.5 median(y(mBs))],[],[],[],[],[-inf -inf  0 -inf],[inf inf 1 inf],[],opt);
        % baseline fit
        bsFit = fitFun(coef,t);
     

    case 'spline'
        n_knots = params(1);
        splOrd = params(2);  % spline order

        spl = spap2(n_knots, splOrd, t(:), y(:), double(mBs(:)));
        % newknt for a possibly better knot distribution
        knots = newknt(spl);
        % least-squares approximation
        spl = spap2(knots, splOrd, t(:), y(:), double(mBs(:)));
        bsFit = fnval(spl, t);

    case 'SmoothingSpline'
        s_span = params(1);
        ft = fittype('smoothingspline');
        opts = fitoptions('Method','SmoothingSpline',...
            'SmoothingParam',s_span,...
            'Normalize','on',...
            'Exclude',~mBs);
        % fit
        [f_bs, ~,~] = fit(t(:), y(:), ft, opts);
        bsFit = feval(f_bs, t);


    otherwise
        ft = fittype(type);
        opts = fitoptions('Method','NonlinearLeastSquares',...
            'Normalize','on',...
            'Robust','Bisquare',...
            'Exclude',~mBs);
        % fit
        [f_bs, ~,~] = fit(t(:), y(:), ft, opts);
        bsFit = feval(f_bs,t);
        
end

if plotResult
    % delete previous fit
    delete(findobj(ax,'Tag','baselineFit'))
    % new fit
    line(t, bsFit, 'Parent',ax, ...
        'Color','b', 'LineStyle','-', 'LineWidth',2,...
        'Tag','baselineFit');
    
    % delete previous mask
    hl_m_Old = findobj(ax,'Type','Line','-regexp','Tag','Mask');
    delete(hl_m_Old)
    % new mask
    hl_m = line(t(mBs),y(mBs), 'Parent',ax, 'Color','g',...
        'LineStyle','none', 'Marker','.', 'MarkerSize',20, ...
        'LineWidth',1, 'Tag','baselineMask');
    uistack(hl_m, 'bottom')
    
end

end

