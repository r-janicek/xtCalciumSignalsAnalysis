function coefEst = paramEstimEMG(y,t,t0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% estimation of parameters of exponentially modified gaussian (EMG)
% y = data
% t = time
% t0 = beggining of EMG
keyboard
[~,ind_t0] = min( abs(t-t0) );

% EMG only
yEMG = y(ind_t0:end);
tEMG = t(ind_t0:end);


mean(yEMG./max(yEMG).*tEMG)

% position of mean
[~,p_m] = min( abs(yEMG - mean(yEMG))); 
p_m = tEMG(p_m);

% position of median
[~,p_med] = min( abs(yEMG - median(yEMG))); 
p_med = tEMG(p_med);


% estimate skewnewss of EMG
s = abs((mean(yEMG) - median(yEMG)))/std(tEMG,yEMG);

% estimation of parameters of EMG
m = mean(yEMG) - std(yEMG)*(s/2)^(1/3)
sd = std(yEMG)^2 * (1-(s/2)^(2/3));
tau = std(yEMG) * (s/2)^(2/3);

figure
plot(tEMG,yEMG)

end

