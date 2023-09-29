function baselineCorrection(~,~,main_fig)

set(getappdata(main_fig,'main_fig'),'Pointer','watch')
drawnow
keyboard
% baseline correction
ImageData = getappdata(main_fig,'Img_data');
x_t_img = getappdata(main_fig,'xt_img'); 
ax_img = getappdata(main_fig,'ax_Img');
ax_prof = getappdata(main_fig,'ax_prof');
pxSz = getappdata(main_fig,'pxSize_x');

n_px = round(1/pxSz); 

if mod(n_px,2)==0    
    n_px = n_px + 1;   
end
d = floor(n_px/2);

% calculate break points for detrending
N_intervals = round(x_t_img(end)/1000);
bp = (1:1:N_intervals-1).*(round(numel(x_t_img)/N_intervals));
bp = [bp,numel(x_t_img)];

% Set up fittype and options.
ft_pb = fittype('poly2'); %parts (1000 ms) of baseline
ft_wb = fittype('poly5'); %whole baseline
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Normalize = 'on';
opts.Robust = 'Bisquare';

t = linspace(1,size(ImageData,2),size(ImageData,2));
Img_c = zeros(size(ImageData));

ImageData_s = ImageData;           
ImageData_s = [ImageData_s(1:d,:);ImageData_s;ImageData_s(end-d+1:end,:)];

for i=d+1:size(ImageData_s,1)-d
    
    % take average from 1 um to estimate trend of baseline
    y = mean(ImageData_s(i-d:i+d,:),1); 
%     y_raw = mean(ImageData_s(i,:),1);
%      keyboard  
%     a=0;
%     for j=1:numel(bp)
%         
%         t_pb = t(1+(j-1)*bp(j-a):bp(j));
%         y_pb = y_avrg(1+(j-1)*bp(j-a):bp(j));
%         
%         [f_p, ~,~] = fit(t_pb(:),y_pb(:), ft_pb, opts );
%         fit_pb = feval(f_p,t_pb(:));
%         a = a+1;
%         
%         figure
%         plot(t_pb,y_pb)
%         hold on 
%         plot(t_pb,fit_pb,'r')
%         
%         
%     end
%     
%     
%     % detrend, linear
%     y_d = detrend(y_avrg,'linear',bp);
%     y_line = y_avrg-y_d;  
%     % subtract trend from original signal
%     y = y_raw - (y_line -  mean( y_raw(y_raw>=1.5*std() )  )           );
%     
    
    
    
    % fit signal with subtracted trend with high order polynom     
    % exclude events > 75th percentile      
    perc = prctile(y,[50 75]);    
    m_exc = y > perc(2);
    
    opts.Exclude = m_exc;
    %opts.Weights = ~m_exc*10;
    % Fit model to data.
    [f, ~,~] = fit(t(:),y(:), ft_wb, opts );
    p_ss = feval(f,t);
    
%   
% %     
%    keyboard
%  
%    tic
%    figure
%    plot(y_raw)
%    hold on 
%    plot(smooth(y_raw,0.2,'rlowess'),'r')
%    toc
   
%      
%     
%     
    
    Img_c(i,:) = ImageData_s(i,:) - (p_ss' - min(p_ss(:)));
               
end

%remove zeros
Img_c(1:d,:)=[];

  
set(findobj(ax_img,'Type','Image'),'CData',Img_c)

set(ax_prof.Children,'XData',x_t_img,'YData',mean(Img_c,1))
set(ax_prof,'YLim',[min(mean(Img_c,1)) max(mean(Img_c,1))])

cla(getappdata(main_fig,'ax_img_sparks'));
cla(getappdata(main_fig,'ax_img_sparks_2'));

setappdata(main_fig,'Img_data',Img_c);

set(getappdata(main_fig,'main_fig'),'Pointer','arrow')
drawnow

end

