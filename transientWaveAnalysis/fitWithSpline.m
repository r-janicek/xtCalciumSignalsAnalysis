function [splFit,splFitUps,t_ups] = fitWithSpline(t,y,nK,splOrd,pxSzT,coefFitEventRise)

pxT_ups = 0.1;   % ms
t_ups = (t(1):pxT_ups:t(end));

% check if there are parameters from exp rise
if isempty(coefFitEventRise)
    
    % first guess of spline coefficients
    spl0  = spap2(nK, splOrd, t, y);
    % newknt for a possibly better knot distribution
    knots = newknt(spl0);
    spl0  = spap2(knots, splOrd, t, y);
    
    splFit = fnval(spmak(spl0.knots,spl0.coefs),t);
    splFitUps = fnval(spmak(spl0.knots,spl0.coefs),t_ups);
    
else
    % params from event rise fit
    bs = coefFitEventRise(:,4);
    t0 = coefFitEventRise(:,1);
    
    % fit parts of profile with spline
    nKp = ceil(nK/numel(t0));
    
    % separate parts with individuals peaks 
    splFit = zeros(size(t));
    splFitUps = zeros(size(t_ups));
    for i = 1:numel(t0)
           
        if i==numel(t0)
            % start of event in points (in t_ups)
            [~,eS] = min( abs( t_ups - t0(i) ) );
            eE = numel(t_ups);
            % in real time axis, in points
            [~,eSr] = min( abs( t - t0(i) ) );
            eEr = numel(t);
            bsP = 0; % baseline
            
        else
            % start of event in points (in t_ups)
            [~,eS] = min( abs( t_ups - t0(i) ) );
            [~,eE] = min( abs( t_ups - t0(i+1) ) );
            eE = eE - 1;
            % in real time axis, in points
            [~,eSr] = min( abs( t - t0(i) ) );
            [~,eEr] = min( abs( t - t0(i+1) ) );
            eEr = eEr - 1;
            bsP = 0; % baseline
            
        end
        
        yP = y(eSr:eEr);
        tP = t(eSr:eEr);
        tupsP = t_ups(eS:eE);

        % first guess of spline coefficients
        spl0  = spap2(nKp, splOrd, tP, yP);
        % newknt for a possibly better knot distribution
        knots = newknt(spl0);
        spl0  = spap2(knots, splOrd, tP, yP);
        
        splFitP = fnval(spmak(spl0.knots,spl0.coefs),tP);
        splFitUpsP = fnval(spmak(spl0.knots,spl0.coefs),tupsP);
        
        splFit(eSr:eEr) = splFitP;
        splFitUps(eS:eE) = splFitUpsP;
        
        if i==1
            splFit(1:eSr) = bs(1);
            splFitUps(1:eS) = bs(1);
        end
        
    end   
     
end

