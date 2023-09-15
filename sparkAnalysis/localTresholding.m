function [ output_args ] = localTresholding( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%%%%%%%%%%%%%
% method from http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5695425
keyboard

I = crop(:,1:300);
I = imgaussfilt(I,1);

std_img = std2(I);
tresh_std = 0.3*std_img;

% find radius for each point x,y
r_max = 10;
[xx, yy] = meshgrid(1:size(I,2),1:size(I,1));
I_n = zeros(size(I));
tic
for y = 1:size(I,1)
     
    for x = 1:size(I,2)
       
        % get part of image        
        cc = x-r_max:x+r_max;
        cc(cc<=0)=[];
        cc(cc>size(I,2))=[];
        rr = y-r_max:y+r_max; 
        rr(rr<=0)=[];
        rr(rr>size(I,1))=[];
        
        img_p = I(rr,cc);
        
        cond=1;       
        r=1;
        while cond && r<=r_max
            
            C = sqrt((xx(rr,cc)-x).^2+(yy(rr,cc)-y).^2)<=r;           
            l_std = std(img_p(C));
            
            if l_std > tresh_std
                cond = 0;                
            else
                r=r+1;
            end
        end
        
        l_mean = mean(img_p(C));
        l_std = max(l_std,tresh_std);
        
        I_n(y,x) = (I(y,x) - l_mean) / l_std;
       
        %clearvars img_p 
    end
    
end
toc

figure
sp1=subplot(2,1,1);
image(I,'CDataMapping','scaled')
sp2=subplot(2,1,2);
image(I_n,'CDataMapping','scaled')

linkaxes([sp1,sp2],'xy')


%%%%%%%%%%%%%


end

