function bsFit = fitBaseline(t,y,type,mBs,plotResult,ax)
%keyboard
% t = time
% y = profile
% type = type of function
% mBs = mask of baseline
% plotResult = 1 or 0, plot results in axes ax

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
        bsFit = feval(f_bs,t);
        
            
    case 'stretchedExp'
        fitFun = @(p,x) p(1).*exp( -(x./p(2)).^p(3) ) + p(4);
        fitFunSum = @(p,t,y) sum( (y-fitFun(p,t)).^2 );
        % do fit
        opt = optimoptions('fmincon','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,...
            'MaxIter',1000,'MaxFunEvals',3000);
        coef = fmincon(@(p)fitFunSum(p,t(mBs),y(mBs)),...
            [100 t(end)/10 0.5 median(y(mBs))],[],[],[],[],[-inf -inf  0 -inf],[inf inf 1 inf],[],opt);
%         mean(prof.y(find(prof.baselineM,100,'last')))
%            figure
%            plot(prof.t,prof.y)
%            hold on
%            plot(prof.t,fitFun(coef,prof.t),'r')
%            plot(prof.t,fitFun([1 100 0.25 median(prof.y(prof.baselineM))],prof.t),'g')
                
        % baseline fit
        bsFit = fitFun(coef,t);
        
        
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
    line(t,bsFit,'Parent',ax,'Color','b','LineStyle','-','LineWidth',2,...
        'Tag','baselineFit');
    
    % delete previous mask
    hl_m_Old = findobj(ax,'Type','Line','-regexp','Tag','Mask');
    delete(hl_m_Old)
    % new mask
    hl_m = line(t(mBs),y(mBs),'Parent',ax,'Color','g',...
        'LineStyle','none','Marker','.','MarkerSize',20,'LineWidth',1,'Tag','baselineMask');
    uistack(hl_m, 'bottom')
    
end

end

