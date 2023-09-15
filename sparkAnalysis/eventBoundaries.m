function out = eventBoundaries(data,eventsPos,main_fig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% data = image data or profile 
% eventPos = image data = structure of detected events 
%          = profile = positions of peaks
keyboard
N_events = numel(eventsPos);

pxSz_t = getappdata(main_fig,'pxSize_t'); % ms
pxSz_x = getappdata(main_fig,'pxSize_x'); % um

% max duration of baseline 
maxDurOfBaselineT = round(100/pxSz_t); % in points
maxDurOfBaselineX = round(5/pxSz_t); % in points

[XOut, YOut, ZOut] = prepareSurfaceData((1:1:size(data,2)), (1:1:size(data,1)), data);
f=fit([XOut, YOut], ZOut,'loess');

fit_data = csaps({1:1:size(data,1),1:1:size(data,2)},data,val,{1:1:size(data,1),1:1:size(data,2)});

% check dinemsions of data
Ndim = ndims(data);

switch Ndim
    case 1 % profile
        

    case 2 % image
        
             
end

csaps([XOut, YOut], ZOut)


startOfSpark = zeros(length(pks),1); 
endOfSpark = zeros(length(pks),1);
detectedEventsMask = false(numel(prof_t),1);




% smooth profile for difference analysis, loess with defined duration in ms  
%prof_s = smooth(prof_t,3);
n_pts = round(smooth_span/pxSz_t);
prof_s = smooth(prof_t,n_pts/length(prof_t),'loess');

bs_crit = round(bs_crit);

%remove baseline, only events bigger than specified percentile stay
percentl = prctile(prof_s,[25 50 bs_crit]);
%iqr = percentl(3)- percentl(1);

prof_s(prof_s < percentl(3))=percentl(3);

% figure
% plot(prof_s)


[v_s,l_s] = findpeaks(prof_s);

    
%     % do not do extension for sparks detection
%     if ~strcmp(analysis_type,'spark detection')
%         
%         ext_x = round(2/pxSz);   % extension in x, in pixels        
%         ext_t_s = cols_e(1) - sparkStart(i); % extension in t, in pixels
%         
%         % set the extension of end of spark in pixels, fixed to 100 ms if more
%         % or last one
%         if i==length(S_event)
%             ext_t_e = round(100/pxSize_t);
%         else
%             ext_t_e = sparkStart(i+1)-cols_e(end);
%             
%             if ext_t_e*pxSize_t > 100
%                 
%                 ext_t_e = round(100/pxSize_t);
%                 
%             end
%         end
%         
%         %check if it is OK in terms of size of Img_data
%         if (size(Img_data,1)<rows_e(end)+ext_x) || (rows_e(1)-ext_x <= 0)
%             
%             ext_x = min(size(Img_data,1)-rows_e(end)-1, rows_e(1)-1);
%             
%         end
%         
%         if (cols_e(1)-ext_t_s) <= 0
%             
%             ext_t_s = cols_e(1)-1;
%             
%         end
%         
%         if size(Img_data,2) < (cols_e(end)+ext_t_e)
%             
%             ext_t_e = size(Img_data,2)-cols_e(end)-1;
%             
%         end
%         
%         %expand spark area
%         rows_e = [(rows_e(1)-ext_x:rows_e(1)-1),rows_e,(rows_e(end)+1:rows_e(end)+ext_x)];
%         cols_e = [(cols_e(1)-ext_t_s:cols_e(1)-1),cols_e,(cols_e(end)+1:cols_e(end)+ext_t_e)];
%     end




keyboard




end

