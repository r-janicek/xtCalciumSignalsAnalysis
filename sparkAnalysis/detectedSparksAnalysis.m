function detectedSparksAnalysis(h_O,~,mainFig)

% get data
hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');
sparkDetection = getappdata(mainFig,'sparkDetection');

% selected method of parameters calculation
switch hObjs.calcSpParamMethodRBgroup.SelectedObject.String
    
    case '2D gaussian fit'
        calcMethod = '2DGauss';
        
    case 'max crossing profiles'
        calcMethod = 'peakXTProfiles';
        
    case 'estimate from event img'
        calcMethod = 'estimate from event img';
        
end

% first check if sparks were detected
if isempty(sparkDetection)
    if ~isfield(sparkDetection,'detectedEventsRec')
        errordlg('FIRST DETECT SPARKS!')
        return
    end
end

% remove old analysis
if isfield(sparkDetection,'eventParams')
    sparkDetection = rmfield(sparkDetection,'eventParams');
end
if isfield(sparkDetection,'maskOfAcceptedSparks')
    sparkDetection = rmfield(sparkDetection,'maskOfAcceptedSparks');
end


% get sparks parameters
sparkDetection.eventParams = findDetectedSparksParams(imgData.imgDataXTfluoR, ...
    sparkDetection.detectedEvents, mainFig, calcMethod, sparkDetection.detectedEventsRec);

% save data
setappdata(mainFig,'sparkDetection',sparkDetection);

% sparks filtering based on their parameters, also change colors of rectangles of rejected
% sparks
sparkParamsFiltering([],[],mainFig);


% set(getappdata(mainFig,'main_fig'),'Pointer','watch')
% drawnow
% 
% % dat na vyber analyzu v popup menu a podla toho vyskladat layout 
% if ~getappdata(mainFig,'norm_flag')
%     warndlg('normalize data');
%     return
% end
% 
% % get data
% Img_data = double(getappdata(mainFig,'Img_data'));
% pxSize_t = getappdata(mainFig,'pxSize_t');
% 
% %get type of analysis
% analysis_type = getappdata(mainFig,'analysisType');
%  
% 
% if strcmp(analysis_type,'spark recovery') || ...
%    (strcmp(analysis_type,'spark detection') && strcmp(get(h_O,'Tag'),'findEvents'))
%         
%     % clear previous
%     if isfield(getappdata(mainFig),'S_event')
%         
%         if isfield(getappdata(mainFig),'spark_rec')
%             delete(getappdata(mainFig,'spark_rec'))
%             rmappdata(getappdata(mainFig,'main_fig'),'spark_rec')
%         end
%         
%         if isfield(getappdata(mainFig),'mask_massRec')
%             rmappdata(getappdata(mainFig,'main_fig'),'mask_massRec')
%             rmappdata(getappdata(mainFig,'main_fig'),'crit_mass')
%         end
%         
%         rmappdata(getappdata(mainFig,'main_fig'),'S_event');
%         rmappdata(getappdata(mainFig,'main_fig'),'S_event_recAll');
%         
%     end
%     
%     S_event = findEvents(mainFig,Img_data);
%     %plot bounding rectangles
%     Img_data_events1 = zeros(size(Img_data,1),size(Img_data,2));
%     for i = 1:length(S_event)
%         
%         %Img_data_events1(S_event_filt(i).SubarrayIdx{1},S_event_filt(i).SubarrayIdx{2}) = Img_data_raw(S_event_filt(i).SubarrayIdx{1},S_event_filt(i).SubarrayIdx{2});
%         Img_data_events1(S_event(i).PixelIdxList) = Img_data(S_event(i).PixelIdxList);
%         pos = S_event(i).BoundingBox;
%         spark_rec(i,1) = rectangle('Position',[pos(1)*pxSize_t pos(2) pos(3)*pxSize_t pos(4)],...
%             'Parent',getappdata(mainFig,'ax_Img'),'EdgeColor','r');
%         
%         clearvars pos
%     end
%     
%     setappdata(mainFig,'S_event',S_event);
%     setappdata(mainFig,'S_event_recAll',S_event);
%     setappdata(mainFig,'spark_rec',spark_rec);
%     setappdata(mainFig,'lastPressedPusbutton',h_O);
%       
%     % calculate frequency of sparks
%     calcSparkFreq(mainFig)
%     
% else    
%     if ~isfield(getappdata(mainFig),'S_event')
%         set(getappdata(mainFig,'main_fig'),'Pointer','arrow')
%         drawnow
%         return
%     end
%  
%     S_event = getappdata(mainFig,'S_event');  
%     eventParams = calcEventParams(Img_data,S_event,[],[],[],mainFig);               
%     setappdata(mainFig,'calciumSparksParameters',eventParams);
%        
% end
% 
% 
% % longLasting = get(getappdata(main_fig,'check_longLasting'),'Value');
% 
% set(getappdata(mainFig,'main_fig'),'Pointer','arrow')
% drawnow

end



