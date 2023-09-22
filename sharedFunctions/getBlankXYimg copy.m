function getBlankXYimg(~, ~,mainFig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

imgData = getappdata(mainFig,'imgData');
imgDataXY = imgData.imgDataXY;

sc = get(0,'ScreenSize');
r = sc(3)/sc(4);

hf = figure('Name','get blank from xy image','units','normalized','outerposition',[0.2 0.2 0.4 0.4*r]);
set(hf, 'PaperPositionMode', 'auto','PaperOrientation', 'landscape','PaperType', 'A4');

ha = axes('Parent',hf,'Position',[0.15 0.25 0.7 0.7]); 
image(imgDataXY{1,1},'YData',[1 size(imgDataXY{1,1},1)],'XData',[1 size(imgDataXY{1,1},2)],...
                'CDataMapping','scaled','Parent',ha);
colormap('jet')

scanLinePos = imgData.scanLinePos;
TPPpointPos = imgData.TPPpointPos;
d_TP_laser = imgData.d_TP_laser;
pxSzX = imgData.pxSzX;

if ~isempty(scanLinePos)
    line([scanLinePos(1) scanLinePos(2)],[scanLinePos(3) scanLinePos(4)],'Parent',ha,'LineWidth',2,'Color','y',...
        'LineStyle','-');
end

if ~isempty(TPPpointPos)
    pos = [TPPpointPos(1)-(d_TP_laser/pxSzX)/2 TPPpointPos(2)-(d_TP_laser/pxSzX)/2 (d_TP_laser/pxSzX) (d_TP_laser/pxSzX)];
    rectangle('Position',pos,'Curvature',[1 1],'Parent',ha,'EdgeColor','r','LineWidth',2)
end

% scale 10 um
line([ha.XLim(2)-5-10/pxSzX ha.XLim(2)-5],[ha.YLim(2)-5 ha.YLim(2)-5],'Parent',ha,'LineWidth',2,'Color','g',...
    'LineStyle','-');
text(ha.XLim(2)-5-5/pxSzX,ha.YLim(2)-10,'10 \mum','HorizontalAlignment','center',...
    'Color','g','FontSize',12,'Parent',ha,'FontWeight','bold')


h_rect = imrect(ha,[5 5 size(imgDataXY{1,1},2)/4 size(imgDataXY{1,1},1)/4]);
setColor(h_rect,'r')

% constrains for ROI
fcn = makeConstrainToRectFcn('imrect',get(ha,'XLim'),get(ha,'YLim'));
setPositionConstraintFcn(h_rect,fcn); 
setResizable(h_rect,1);

chNames = arrayfun(@(x) sprintf('ch #%d',x) ,(1:numel(imgDataXY(:,1))), 'UniformOutput',0);

uicontrol('Style', 'pushbutton','String','get blank',...
    'FontUnits','normalized','Parent',hf,'Units','normalized','FontSize',0.2,...
    'Position', [0.15 0.05 0.15 0.15],'Callback', {@getBlankYX,mainFig,h_rect,hf});
   
uicontrol('Style','popup',...
    'Parent',hf,'Units','normalized','Position', [0.4 0.05 0.15 0.1],...
    'FontUnits','normalized','FontSize',0.2,'Tag','xyImg','Value',1,...
    'String',chNames,'Callback', {@setImgChannel,mainFig,ha});

end

