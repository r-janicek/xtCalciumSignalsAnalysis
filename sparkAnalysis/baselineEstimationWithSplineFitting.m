function [ output_args ] = baselineEstimationWithSplineFitting( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

keyboard
im = double(getappdata(main_fig,'rawImageData'));

% (im, filt_space_fwhm, smooth_iter, base_iter, do_vst, tau, expansion_factor)
% (im, 8, 3, 5, 1, 5, 0.5)


[bw, im_detect] = spark_detect_vst(im, 8, 3, 5, 1, 5, 0.25);

figure
subplot(2,1,1)
imagesc(crop)
subplot(2,1,2)
imagesc(bw);


pxSize_t = getappdata(main_fig,'pxSize_t');
pxSz = getappdata(main_fig,'pxSize_x');
x_t = getappdata(main_fig,'xt_img');
ax_img_sparks = getappdata(main_fig,'ax_img_sparks');
ax_img_sparks_2 = getappdata(main_fig,'ax_img_sparks_2');
crit = str2double(get(getappdata(main_fig,'h_edit_tresh_l'),'String'));
smooth_span = str2double(get(getappdata(main_fig,'h_smooth_edit'),'String'));
bs_crit = str2double(get(getappdata(main_fig,'h_bsDet_edit'),'String'));



F = mean(im,1);

% smooth profile with loess with defined duration in ms  
n_pts = round(50/pxSize_t);

p_s = smooth(F,n_pts/length(F),'loess');

t = (1:1:size(im,2));

%[X,Y] = ndgrid(1:1:19,2:2:12)


% G = csaps(t,F,[],t,[]);
% % weights
% w = zeros(size(G));
% w((F - G)<0) = 100;
% w((F - G)>0) = 0.1;
% H = csaps(t,G,[],t,w);


%%
n_knots = ceil(x_t(end)/1000);
splOrd = 3; % spline order, poly3 


x_t = selectedProf.t;
F = selectedProf.y;

spl = spap2(n_knots,splOrd,x_t,F); 
% newknt for a possibly better knot distribution
knots = newknt(spl); 
% least-squares approximation
spl = spap2(knots, splOrd, x_t, F); 
G = fnval(spl,x_t);

% weights
w = zeros(size(G));
w((F - G)<0) = 100;
w((F - G)>0) = 0;
spl = spap2(n_knots,splOrd,x_t,G); 
% newknt for a possibly better knot distribution
knots = newknt(spl); 
% least-squares approximation
spl = spap2(knots, splOrd, x_t, G); 
H = fnval(spl,x_t);

sResInit = sum(abs(G-H));

sRes1 = sResInit;
sRes2 = 0;
i = 0;
while abs(sRes1-sRes2) > 0.01*sResInit
    
    i = i +1
    
    sRes1 = sum(abs(G-H));
    
    Fs = H;
    
    spl = spap2(n_knots,splOrd,x_t,F); 
    % newknt for a possibly better knot distribution
    knots = newknt(spl); 
    % least-squares approximation
    spl = spap2(knots, splOrd, x_t, Fs); 
    G = fnval(spl,x_t);
    
    w((G - Fs)<0) = 1000;
    w((G - Fs)>0) = 0;
    spl = spap2(n_knots,splOrd,x_t,G); 
    % newknt for a possibly better knot distribution
    knots = newknt(spl); 
    % least-squares approximation
    spl = spap2(knots, splOrd, x_t, G); 
    H = fnval(spl,x_t);
       
    sRes2 = sum(abs(G-H));
           
end

figure
plot(F)
hold on
plot(H,'r')


end

