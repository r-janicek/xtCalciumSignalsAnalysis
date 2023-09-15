function selectFitFcn(hO,~,analysisFig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

hObjsA = getappdata(analysisFig,'hObjsA');
popUpFitFun = hObjsA.popUpMenuEventsFcn; 

if ~isempty(hO)
    if ~strcmp(hO.Style,'popupmenu')
        switch hObjsA.check_pieceWise.Value
            
            case 0
                if str2double(hObjsA.h_edit_Npeaks.String)>1
                    set(popUpFitFun,'String',...
                        {'none', 'expModGauss', 'doubleBoltzmann', 'Asym2Sigmoid'});
                    popUpFitFun.Value = 1;
                    
                else
                    set(popUpFitFun,'String',...
                        {'none', 'spline', 'expModGauss', 'doubleBoltzmann', 'Asym2Sigmoid'});
                    popUpFitFun.Value = 1;
                end
                
                
            case 1
                if str2double(hObjsA.h_edit_Npeaks.String)>1
                    set(popUpFitFun,'String',...
                        {'none', '1expR1expD', '1expR2expD', 'CaSpikeFun'});
                    popUpFitFun.Value = 1;
                else
                    set(popUpFitFun,'String',...
                        {'none', 'spline', '1expR1expD', '1expR2expD', 'CaSpikeFun'});
                    popUpFitFun.Value = 1;
                end
                
        end
    end
end

switch popUpFitFun.String{popUpFitFun.Value}
    
    case 'spline'       
        set(hObjsA.h_txt_paramFit1,'String',{'number of'; 'knots'})
        set(hObjsA.h_edit_paramFit1,'String',num2str(ceil(hObjsA.ax_orgImg.XLim(2)/100)),'Enable','on')
        set(hObjsA.h_txt_paramFit2,'String',{'order'; '<2,10>'})
        set(hObjsA.h_edit_paramFit2,'String',num2str(3),'Enable','on')
                
    case 'expModGauss'
        set(hObjsA.h_txt_paramFit1,'String',{'fit param'; '#1:'})
        set(hObjsA.h_edit_paramFit1,'String','nan','Enable','off')
        set(hObjsA.h_txt_paramFit2,'String',{'fit param'; '#2:'})
        set(hObjsA.h_edit_paramFit2,'String','nan','Enable','off')
        set(hObjsA.check_pieceWise,'Value',0)
        
    case 'doubleBoltzmann'
        set(hObjsA.h_txt_paramFit1,'String',{'fit param'; '#1:'})
        set(hObjsA.h_edit_paramFit1,'String','nan','Enable','off')
        set(hObjsA.h_txt_paramFit2,'String',{'fit param'; '#2:'})
        set(hObjsA.h_edit_paramFit2,'String','nan','Enable','off')
        set(hObjsA.check_pieceWise,'Value',0)
    
    case 'Asym2Sigmoid'
        set(hObjsA.h_txt_paramFit1,'String',{'fit param'; '#1:'})
        set(hObjsA.h_edit_paramFit1,'String','nan','Enable','off')
        set(hObjsA.h_txt_paramFit2,'String',{'fit param'; '#2:'})
        set(hObjsA.h_edit_paramFit2,'String','nan','Enable','off')
        set(hObjsA.check_pieceWise,'Value',0)
        
    case '1expR1expD'
        set(hObjsA.h_txt_paramFit1,'String',{'fit param'; '#1:'})
        set(hObjsA.h_edit_paramFit1,'String','nan','Enable','off')
        set(hObjsA.h_txt_paramFit2,'String',{'fit param'; '#2:'})
        set(hObjsA.h_edit_paramFit2,'String','nan','Enable','off')
        set(hObjsA.check_pieceWise,'Value',1)
        
    case 'CaSpikeFun'
        set(hObjsA.h_txt_paramFit1,'String',{'fit param'; '#1:'})
        set(hObjsA.h_edit_paramFit1,'String','nan','Enable','off')
        set(hObjsA.h_txt_paramFit2,'String',{'fit param'; '#2:'})
        set(hObjsA.h_edit_paramFit2,'String','nan','Enable','off')
        set(hObjsA.check_pieceWise,'Value',1)

end


end

