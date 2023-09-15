function calcWholeImgWaveSpeed(~,~,mainFig)


keyboard

%not really working yet

% get data
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');
ax_img = hObjs.ax_img;
pxSzX = imgData.pxSzX;
pxSzT = imgData.pxSzT;

% wave speed estimation 
I = imgData.imgDataXTfluoFN;

figure
imagesc(I)
 
 
% treshold image
% baseline value
mI = false(size(I));
for i=1:size(I,1)
    mL = I(i,:)<(mean(I(i,:))+1*std(I(i,:)));
    bs = mean( I(i,mL) );
    crit = 1.5*std( I(i,mL) );
    mI(i,:) = I(i,:)>(bs+crit);
end

% remove objects with pixels less than: 10um * 100 ms 
nPix = round( (10/pxSzX) * (100/pxSzT) );
BW = bwareaopen(mI,nPix,8);
% fill holes
BW = imfill(BW,'holes');

I(~BW) = 0;

figure
imagesc(I)

% normalized cross correlation
t_I = [0 size(I,2)*pxSzT];
x_I = [0 size(I,1)*pxSzX];

crr = normxcorr2(I,I);
t_crr = [-size(I,2)*pxSzT size(I,2)*pxSzT];
x_crr = [-size(I,1)*pxSzX size(I,1)*pxSzX];

% remove lowest 10 percent
mCrr = crr>0.1;

figure
imagesc(mCrr)

% find all components 
CC = bwconncomp(mCrr);
stats = regionprops(CC,crr,'MaxIntensity','BoundingBox','PixelList','PixelIdxList','SubarrayIdx','Orientation');

fitAreaStat = stats( find( [stats.MaxIntensity] == max( [stats.MaxIntensity] ) ) );

mm = false(size(crr));
mm(fitAreaStat.PixelIdxList)=true;

% get only maximum
crrMax = crr.*mm;

% normalized cross correlation
t_crrMax = (-(size(I,2)-1)*pxSzT : pxSzT : (size(I,2)-1)*pxSzT);
x_crrMax = (-(size(I,1)-1)*pxSzX : pxSzX : (size(I,1)-1)*pxSzX);

t_crrMax = t_crrMax(fitAreaStat.SubarrayIdx{2});
x_crrMax = x_crrMax(fitAreaStat.SubarrayIdx{1});

% 
crrMax = crrMax(fitAreaStat.SubarrayIdx{1},fitAreaStat.SubarrayIdx{2});



figure
imagesc(crrMax)



dt = (t_crrMax(end)-t_crrMax(1))/1000; % in s
dx = x_crrMax(end)-x_crrMax(1);

waveSpeed = dx/dt; % in um/s




s = svd(crrMax);

s(1)





end
