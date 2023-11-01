function checkboxFcn(hO,~,mainFig)

hObjs = getappdata(mainFig,'hObjs');
imgData = getappdata(mainFig,'imgData');

profLineObj = findall(hObjs.ax_prof, 'Type','Line', ...
    'Tag','wholeCellProfile');
prof = profLineObj.YData;
t = profLineObj.XData;

switch hO.Tag
    
    case 'check_bsFit'
        switch hO.Value
            case 0
                hObjs.popUpMenuBaselineFcn.Enable = 'off';
                % delete previous baseline fit line
                delete(findobj(hObjs.ax_prof, 'Type','Line', ...
                    'Tag','baselineFit'))
                % delete previous baseline mask line
                delete(findobj(hObjs.ax_prof, 'Type','Line', ...
                    '-regexp','Tag','Mask'))
            case 1
                hObjs.popUpMenuBaselineFcn.Enable = 'on';
                % check if there is an old baseline mask
                if isfield(imgData, 'baselineM')
                    mBs = imgData.baselineM;
                    if sum(imgData.baselineM)<10
                        mBs = prof <= prctile( prof, hObjs.sld_Bs.Value);
                    end
                else
                    % calculate first mask of baseline
                    mBs = prof <= prctile( prof, hObjs.sld_Bs.Value);
                end
                if strcmp(hObjs.popUpMenuBaselineFcn.String{...
                        hObjs.popUpMenuBaselineFcn.Value},'spline')
                    mBs(1) = true;
                    mBs(end) = true;
                end
                
                imgData.baselineM = mBs(:);
                setappdata(mainFig, 'imgData', imgData);
                
                % do fitting of baseline
                type = hObjs.popUpMenuBaselineFcn.String{...
                    hObjs.popUpMenuBaselineFcn.Value};
                bsFitOut = fitBaseline(t, prof, ...
                    type, mBs, 1, hObjs.ax_prof);
                imgData.bsFit = bsFitOut;
                setappdata(mainFig, 'imgData', imgData);                 
        end
        
    otherwise
        return
               
end

