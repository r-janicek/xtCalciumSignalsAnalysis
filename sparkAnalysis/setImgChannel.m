function setImgChannel(hO,~,mainFig,xyImgAx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

hObjs = getappdata(mainFig,'hObjs');

set(mainFig,'Pointer','watch')
drawnow

popUpVal = hO.Value;
popUpStr = hO.String;
ch = popUpStr{popUpVal};
ch_num = regexp(ch,'\d*','match');
try
    ch_num = str2double(ch_num{1});
catch
    ch_num = [];
end
imgData = getappdata(mainFig,'imgData');

switch hO.Tag
      
    case 'mainImgFluoCh'
           
        if isempty(ch_num)            
            errordlg('Fluo channel cannot be empty!');  
            hO.Value = 1;
            return
        end
        
        pxSzT = imgData.pxSzT;
        pxSzX = imgData.pxSzX;
        
        imgData.wholeImgFluoXT = imgFiltering(imgData.wholeImgXT{ch_num,1},pxSzT,pxSzX);
        imgData.imgDataXTfluoR = imgData.imgDataXT{ch_num,1};
        imgData.imgDataXTfluoRN = imgData.imgDataXT{ch_num,1};
        imgData.imgDataXTfluoF = imgFiltering(imgData.imgDataXT{ch_num,1},pxSzT,pxSzX);
        imgData.imgDataXTfluoFN = imgData.imgDataXTfluoF;
        
        % look if now trans channel is the same number 
        if hObjs.popUpMenuChTrans.Value == popUpVal
            
            numChs = regexp(popUpStr,'\d*','match');
            numChs = cellfun(@(x) str2double(x),numChs,'UniformOutput',0);
            
            m1 = cellfun(@(x) ~isempty(x),numChs,'UniformOutput',1);
            m2 = cellfun(@(x) ~isequal(x,popUpVal),numChs,'UniformOutput',1); 
            
            transChVal = find((m1&m2)==1,1,'first');
            
            % set new value of trans popupmenu
            hObjs.popUpMenuChTrans.Value = transChVal;
            
            % save new trans channel
            if transChVal > numel(imgData.imgDataXT)
                imgData.imgDataXTtrans = [];             
            else
                imgData.imgDataXTtrans = imgData.imgDataXT{transChVal,1};
            end
            
        end
      
        % show images of newly assigned channels
        h_imgTrans = findobj(hObjs.h_ax_transCh,'Type','image');  
        h_imgTrans.CData = imgData.imgDataXTtrans;
        set(hObjs.h_ax_transCh,'XLim',[0.5 size(imgData.imgDataXTtrans,2)+0.5],...
                               'YLim',[0.5 size(imgData.imgDataXTtrans,1)+0.5])
                            
        h_imgFluo = findobj(hObjs.ax_img,'Type','image');
        h_imgFluo.CData = imgData.imgDataXTfluoF;
        set(hObjs.ax_img,'XLim',[imgData.t(1) imgData.t(end)],...
                         'YLim',[0.5 size(imgData.imgDataXTfluoF,1)+0.5])
      
        % time profile
        prof_t = mean(imgData.imgDataXTfluoF,1);             
        h_l = get(hObjs.ax_prof,'Children');
        set(h_l,'YData',prof_t);
        try
            ylim(hObjs.ax_prof,[min(prof_t)*0.95 max(prof_t)*1.05])
        catch
            ylim(hObjs.ax_prof,[min(prof_t)-1 max(prof_t)+1])
        end
        setappdata(mainFig,'imgData',imgData);
                       
    case 'mainImgTransCh'
        
        % save new trans channel
        if isempty(ch_num)
            imgData.imgDataXTtrans = [];
        else
            imgData.imgDataXTtrans = imgData.imgDataXT{ch_num,1};
        end  
        
        % show images of newly assigned channels
        h_imgTrans = findobj(hObjs.h_ax_transCh,'Type','image');  
        h_imgTrans.CData = imgData.imgDataXTtrans;
        try
            set(hObjs.h_ax_transCh,'XLim',[0.5 size(imgData.imgDataXTtrans,2)+0.5],...
                'YLim',[0.5 size(imgData.imgDataXTtrans,1)+0.5])
         
            if strcmp(getappdata(mainFig,'analysisType'),'spark recovery photolysis')
                % calculate new positions of photolysis pulses      
                [imgData.s_TPP,imgData.e_TPP,...
                imgData.durOfTPP,imgData.d_TP_laser,...
                imgData.posOfTPPinScanLine] = loadPhotolysisPositions(...
                                                    imgData.imgDataXTfluoF,...
                                                    imgData.imgDataXTtrans,...
                                                    imgData.t,...
                                                    imgData.pxSzT,...
                                                    imgData.pxSzX,...
                                                    imgData.TPPpointPos,...
                                                    imgData.scanLinePos,...
                                                    hObjs.ax_img);
            
                                                
            end
            
        catch
            
        end
        setappdata(mainFig,'imgData',imgData);
        
        
    case 'xyImg'
      
        imgDataXY = imgData.imgDataXY;
        xyImg_ch_d = double(imgDataXY{popUpVal,1});       
        xyImg = findobj(xyImgAx,'Type','image');
        set(xyImg,'CData',xyImg_ch_d);
        
end

set(mainFig,'Pointer','arrow')
drawnow


end

