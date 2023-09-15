function axesTransChFunction(hO,E)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

switch E.Button
    
    case 1 % create window
        
        % check if window exist
        hf = findall(0,'tag','2CH');
        
        if isempty(hf)
                      
            hf = figure('Name','zoom trans channel','Tag','2CH',...
                'units','normalized','outerposition',[0 0.5 1 0.35]);
            ax = axes('Parent',hf);
            set(ax,'Position',[0.03 0.01 0.94 0.98])
            
            switch hO.Type
                
                case 'image'
                    image(hO.CData,'CDataMapping','scaled','Parent',ax);
                    set(ax,'XTick',[],'YTick',[],'PickableParts','all')
                    
                case 'axes'
                    image(hO.Children.CData,'CDataMapping','scaled','Parent',ax);
                    set(ax,'XTick',[],'YTick',[],'PickableParts','all')
            end
            
            
        else
            % Figure exists so bring Figure to the focus
            figure(hf);
        end
        
        
    case 3 %delete window
        
        delete(findobj(0,'Type','figure','Tag','2CH','Name','zoom trans channel'))
        
end


end

