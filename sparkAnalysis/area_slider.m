function area_slider(~,~,main_fig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

typeAnalysis = get(getappdata(main_fig,'lastPressedPusbutton'),'Tag'); 
S_event_SpRec = getappdata(main_fig,'S_event_recAll');
mass = cellfun(@sum ,{S_event_SpRec.PixelValues});


crit_mass = round(get(getappdata(main_fig,'sld_area'),'Value'));

mask_massRec = mass>crit_mass;

if isfield(getappdata(main_fig),'mask_massRec')
    rmappdata(getappdata(main_fig,'main_fig'),'mask_massRec');
end

setappdata(main_fig,'mask_massRec',mask_massRec)
setappdata(main_fig,'crit_mass',crit_mass)

set(getappdata(main_fig,'sld_area'),'Value',crit_mass)
set(getappdata(main_fig,'txt_sld_area'),'String',sprintf('eventPxMass: %g',crit_mass))

switch typeAnalysis  
    
    case 'findEvents'
             
        spark_rec = getappdata(main_fig,'spark_rec');
        arrayfun(@(x) set(x,'Visible','off'),spark_rec);
        arrayfun(@(x) set(x,'Visible','on'),spark_rec(mask_massRec));
         
        setappdata(main_fig,'S_event',S_event_SpRec(mask_massRec));
        calcSparkFreq(main_fig)
        
    otherwise
              
       crosses = [S_event_SpRec.centreLine];                 
       arrayfun(@(x) set(x,'Visible','off'),crosses);         
       arrayfun(@(x) set(x,'Visible','on'),crosses(mask_massRec)); 
         
       setappdata(main_fig,'S_event_SpRec',S_event_SpRec(mask_massRec));
          
%          h_rect_prof_pos = getPosition(getappdata(main_fig,'h_rect_prof'));
%          plotProfile(h_rect_prof_pos,main_fig)
                
end

end

