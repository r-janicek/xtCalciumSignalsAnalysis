function out = compareFits(fit1, fit2, nestedModels)
% fit is structure with sum of squares and number of degrees of freedom
% nestedModels = true or false
% compare fits:
% 1. same model fitted to two different datasets
% 2. two different models to same dataset
% 3. compare parameters of models fitted

% calculate F-ratio for nested models, Akaike criterion 
% and use t-test to compare parameters of models 

% check which model has less params (simpler model = null model)
if fit1.numOfParams < fit2.numOfParams
    nullModel = fit1;
    altModel = fit2;
else
    nullModel = fit2;
    altModel = fit1;
end

if nestedModels
    % calculate F-ratio
    F = ( (nullModel.ss - altModel.ss)/(nullModel.df-altModel.df) ) / ...
        ( altModel.ss/altModel.df );
    % get p value
    out.p_Ftest = fcdf(F, nullModel.df-altModel.df, altModel.df, "upper");
else
    % calculate Akaike's Information Criterion, corrected for low number of
    % observations (N=number of observations, K=number of parameters + 1,
    % SS=sum of squares)
    fun_AICc = @(N,K,SS) N*log(SS/N) +2*K + 2*K*(K+1)/(N-K-1);
    fun_prob_AICc = @(AICc1, AICc2) exp(-0.5*(AICc1-AICc2)) / ...
        (1 + exp(-0.5*(AICc1-AICc2)));
   
    out.AICc_null = fun_AICc(nullModel.numOfObs, nullModel.numOfParams+1, ...
        nullModel.ss);
    out.AICc_alt = fun_AICc(altModel.numOfObs, altModel.numOfParams+1, ...
        altModel.ss);

    out.prob_nullModel = fun_prob_AICc(out.AICc_null, out.AICc_alt);
    out.prob_altModel = fun_prob_AICc(out.AICc_alt, out.AICc_null);
    
    out.evidence_ratio = out.prob_altModel / out.prob_nullModel;

end

end
