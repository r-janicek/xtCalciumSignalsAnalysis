function checkboxFcn(hO,~,mainFig)

hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');

profLineObj = findall(hObjs.ax_prof,'Tag','wholeCellProfile');
prof = profLineObj.YData;
t =  profLineObj.XData;


switch hO.Tag
    
    case 'check_bsFit'
        switch hO.Value
            case 0
                hObjs.popUpMenuBaselineFcn.Enable = 'off';    
            case 1
                hObjs.popUpMenuBaselineFcn.Enable = 'on';     
                % calculate first mask of baseline
                sldVal = hObjs.sld_Bs.Value;
                mBs = prof <= prctile( prof, sldVal );
                
                if strcmp(hObjs.popUpMenuBaselineFcn.String{hObjs.popUpMenuBaselineFcn.Value},'spline')
                    mBs(1) = true;
                    mBs(end) = true;
                end
                
                imgData.baselineM = mBs(:);
                setappdata(mainFig,'imgData',imgData);
                
                % do fitting of baseline
                type = hObjs.popUpMenuBaselineFcn.String{hObjs.popUpMenuBaselineFcn.Value};
                bsFitOut = fitBaseline(t,prof,type,mBs,1,hObjs.ax_prof);
                
                imgData.bsFit = bsFitOut;
                setappdata(mainFig,'imgData',imgData);
         
                
%                 % update baseline slider
% bsSldVal = round(sum(prof.baselineM)/numel(prof.baselineM)*100);
% if bsSldVal==100, bsSldVal=99; end
% if bsSldVal==0, bsSldVal=1; end
% hObjs.h_txt_sld_Bs.String = sprintf('mask of baseline as %d percentile',bsSldVal);
% hObjs.sld_Bs.Value = bsSldVal;
                
                          
        end
        
    otherwise
        return
               
end

