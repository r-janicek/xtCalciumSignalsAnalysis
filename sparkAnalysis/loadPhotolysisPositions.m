function [s_TPP,e_TPP,durOfTPP,TPP_delays,d_TP_laser,posOfTPPinScanLine] = ...
            loadPhotolysisPositions(imgDataXTfluo,imgDataXTtrans,...
                                    t,pxSzT,pxSzX,...
                                    pos_P,pos_l,...
                                    mainImgAx,markPosOfTPP)
    
           
try
    %diameter of TP laser in um
    d_TP_laser = 1.5;
    
    if ~isempty(imgDataXTtrans)
       
        % get position and duration of TP pulse
        avrg_TPP = mean(imgDataXTtrans,1);

        mask_TPP = (avrg_TPP>((max(avrg_TPP)-mean(avrg_TPP))*0.2)+mean(avrg_TPP));
        mask_TPP_dif = diff(mask_TPP);
        
        %start and end TTP in pixels
        s_TPP = (find(mask_TPP_dif>0));
        e_TPP = (find(mask_TPP_dif<0));
        durOfTPP = floor(mean(abs(s_TPP-e_TPP).*pxSzT)/1)*1; %in ms
        
        if durOfTPP == 0
            durOfTPP = 1;
        end
        
        % we used duration of TTPs [1,3,5,10,20,30,50,100...]
        if durOfTPP <= 2
            durOfTPP = 1;
        end
        
        if numel(s_TPP)~=numel(e_TPP)            
            durOfTPP = [];
            s_TPP = [];
            e_TPP = [];
        end
        
    else       
        durOfTPP = [];
        s_TPP = [];
        e_TPP = [];
        
    end
       
    if ~isempty(pos_P) && ~isempty(s_TPP)
        %position of TPP point in line scan, distance in pixels
        %from the start points of scanning line
        try
            % find max in x profile during pulses and use it
            % to determine which position of TPP point in scanning line to take
            [~,m_pos] = max(mean(imgDataXTfluo(:,mask_TPP),2));
            
            pointPosInScanLine = [sqrt((pos_P(1) - pos_l(1))^2 + (pos_P(2) - pos_l(3))^2),...
                sqrt((pos_P(1) - pos_l(2))^2 + (pos_P(2) - pos_l(4))^2)];
            
            [~,pp] = min(abs(m_pos-pointPosInScanLine));
            
            posOfTPPinScanLine = pointPosInScanLine(pp);
            
        catch
            posOfTPPinScanLine = [];
        end
        
    else
        posOfTPPinScanLine = [];
        
    end
       
catch   
    posOfTPPinScanLine = [];
    durOfTPP = [];
    s_TPP = [];
    e_TPP = [];
    d_TP_laser = 0;
end


% calculate delays between photolytic pulses rounded to 5 ms
try
    TPP_dif = round((diff(s_TPP).*pxSzT)./5).*5;
    TPP_dif = TPP_dif(:);
     
    B = unique(TPP_dif);
    Ncount = histc(TPP_dif,B);
    Bs = sortrows([Ncount,B],-1);
    
    if size(Bs,1)==1        
        TPP_delays = [Bs(1,2);nan];
    else       
        TPP_delays = [min(Bs(1:2,2));max(Bs(1:2,2))-min(Bs(1:2,2))];
    end
    
catch       
    TPP_delays = [nan;nan];
end


% mark positions of photolytic pulses
if markPosOfTPP
    
    % remove previous lines of photolysis pulse
    delete(findobj(mainImgAx,'Type','Rectangle','Tag','photolysis'))
    delete(findobj(mainImgAx,'Type','line','Tag','photolysis'))
    
    % plot positions of TPP in axes where main image is
    if ~isempty(s_TPP)
        if ~isempty(posOfTPPinScanLine)
            
            % draw new
            for i=1:numel(s_TPP)
                
                pos = [t(s_TPP(i)), posOfTPPinScanLine - (d_TP_laser/pxSzX)/2, durOfTPP, d_TP_laser/pxSzX];
                rectangle('Position',pos,'Parent',mainImgAx,'EdgeColor','k','LineWidth',2,'Tag','photolysis')
                
            end
            
        else
            
            % draw new
            for i=1:numel(s_TPP)
                
                line([t(s_TPP(i)) t(s_TPP(i))],get(mainImgAx,'YLim'),'Parent',mainImgAx,...
                    'LineWidth',2,'Color','k','LineStyle',':','Tag','photolysis');
                line([t(e_TPP(i)) t(e_TPP(i))],get(mainImgAx,'YLim'),'Parent',mainImgAx,...
                    'LineWidth',2,'Color','k','LineStyle',':','Tag','photolysis');
                
            end
        end
    end
end


end

