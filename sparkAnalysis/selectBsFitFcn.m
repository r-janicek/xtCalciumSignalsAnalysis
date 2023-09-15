function selectBsFitFcn(hO,~,fitFig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

hObjsFit = getappdata(fitFig,'hObjsFit');
imgData = getappdata(fitFig,'imgData'); 

switch hO.String{hO.Value}
    
    case 'spline'       
        set(hObjsFit.h_txt_paramFitBs1,'String','number of knots')
        set(hObjsFit.h_edit_paramFitBs1,'String',num2str(ceil(imgData.t(end)/500)),'Enable','on')
        set(hObjsFit.h_txt_paramFitBs2,'String','order <2,10>')
        set(hObjsFit.h_edit_paramFitBs2,'String',num2str(3),'Enable','on')
                
    case 'SmoothingSpline'
        set(hObjsFit.h_txt_paramFitBs1,'String','span (0,1)')
        s_span = numel(imgData.t)/(numel(imgData.t)+ceil(5/(imgData.pxSzT)));
        set(hObjsFit.h_edit_paramFitBs1,'String',s_span,'Enable','on')
        set(hObjsFit.h_txt_paramFitBs2,'String','fit param #2:')
        set(hObjsFit.h_edit_paramFitBs2,'String','nan','Enable','off')
        
    case 'exp1'
        set(hObjsFit.h_txt_paramFitBs1,'String','fit param #1:')
        set(hObjsFit.h_edit_paramFitBs1,'String','nan','Enable','off')
        set(hObjsFit.h_txt_paramFitBs2,'String','fit param #2:')
        set(hObjsFit.h_edit_paramFitBs2,'String','nan','Enable','off')
        
    case 'exp2'
        set(hObjsFit.h_txt_paramFitBs1,'String','fit param #1:')
        set(hObjsFit.h_edit_paramFitBs1,'String','nan','Enable','off')
        set(hObjsFit.h_txt_paramFitBs2,'String','fit param #2:')
        set(hObjsFit.h_edit_paramFitBs2,'String','nan','Enable','off')
        
    case 'stretchedExp'
        set(hObjsFit.h_txt_paramFitBs1,'String','fit param #1:')
        set(hObjsFit.h_edit_paramFitBs1,'String','nan','Enable','off')
        set(hObjsFit.h_txt_paramFitBs2,'String','fit param #2:')
        set(hObjsFit.h_edit_paramFitBs2,'String','nan','Enable','off')
        
    case 'polynomial'
        set(hObjsFit.h_txt_paramFitBs1,'String','order <1,9>')
        set(hObjsFit.h_edit_paramFitBs1,'String',num2str(4),'Enable','on')
        set(hObjsFit.h_txt_paramFitBs2,'String','spline order')
        set(hObjsFit.h_edit_paramFitBs2,'String','nan','Enable','off')
        
end

% do fitting of baseline
fitBaseline([],[],fitFig)

end

