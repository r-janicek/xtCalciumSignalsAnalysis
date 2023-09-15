function [img_filt,filters] = imgFiltering(img_raw,pxSzT,pxSzX)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% filter image data with median and gaussian filter
% use several steps of filtering with increasing size of neighbourhood, use
% weigted median filter, to preserve structures 

% find max window 
n_t = ceil(5/pxSzT);
n_x = ceil(0.5/pxSzX);

n = min(n_t,n_x);
if mod(n,2)==0 || n == 0
   n = n + 1;  
end

%%%
%n=1;
%%%

if isstruct(img_raw)
    
    %img_filt = structfun(@(x) hmf(x,n),img_raw,'UniformOutput',0);
    img_filt = structfun(@(x) medfilt2(x,[n n],'symmetric'),img_raw,'UniformOutput',0);
    %%img_filt = structfun(@(x) imgaussfilt(x,1,'Padding','symmetric'),img_filt,'UniformOutput',0);
       
elseif iscell(img_raw)
 
    %img_filt = cellfun(@(x) hmf(x,n),img_raw,'UniformOutput',0);
    img_filt = cellfun(@(x) medfilt2(x,[n n],'symmetric'),img_raw,'UniformOutput',0);
    
else
        
    %img_filt = hmf(img_raw,n);
    img_filt = medfilt2(img_raw,[n n],'symmetric');

end

filters = {sprintf('median filter, [%d %d]',n,n)};

end

