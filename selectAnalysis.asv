function selectAnalysis(hO,~)

% get current folder
currentFolder = pwd;
%addpath(genpath(currentFolder))
addpath(currentFolder) % do not include subfolders

switch hO.String
        
    case 'calcium sparks detection'
        
        close all   
        addpath(genpath([currentFolder,'/sparkAnalysis'])) % add path
        % start main window
        MainWindowSparksAndProfiles('spark detection');
        hf = findall(0,'tag','sparkAnalysis');
        % save type of analysis       
        setappdata(hf,'analysisType','spark detection'); 
        % set SparkAnalysis window  
        hObjs = getappdata(hf,'hObjs');
                    
        set(hObjs.h_table_notes,...
            'Data',{'ctrl';'permeabilized myocytes';'50 uM fluo-3-5K';'';''})
        
    
    case 'calcium sparks recovery'    
        
        close all             
        addpath(genpath([currentFolder,'/sparkAnalysis'])) % add path
        MainWindowSparksAndProfiles('spark recovery ryanodine');
        
        hf = findall(0,'tag','sparkAnalysis');
    
        % save type of analysis       
        setappdata(hf,'analysisType','spark recovery ryanodine');
                 
        % set SparkAnalysis window  
        hObjs = getappdata(hf,'hObjs');
        set(hObjs.h_edit_ROI,'Enable','on')
        set(hObjs.h_edit_averageWidth,'Enable','on') 
        set(hObjs.h_edit_tresh_prof,'Enable','on') 
        set(hObjs.h_smooth_edit,'Enable','on') 
        set(hObjs.h_table_profs,'Enable','on') 
        set(hObjs.h_bsDet_edit,'Enable','on') 
        set(hObjs.h_pb_profROI,'Enable','on')
        set(hObjs.h_pb_getPosOfProf,'Enable','on')
                              
        set(hObjs.h_pb_sparksAnalysis,'Enable','off') 
        set(hObjs.h_pb_profilesFit,'Enable','on') 
        set(hObjs.h_pb_profileAnalysis,'Enable','on')       
        set(hObjs.check_pairAnalysis,'Enable','off','Value',0)         
        
        set(hObjs.check_sparkParams,'Enable','on') 
        set(hObjs.check_showEventsFigs,'Enable','on') 
        
        set(hObjs.check_saveProfsAndFits,'Enable','on') 
        set(hObjs.check_saveEventsFigs,'Enable','on') 
        
        set(hObjs.h_pb_dataPreview,'Enable','on') 
          
        set(hObjs.h_table_notes,...
            'Data',{'ctrl';'5 uM fluo-3 AM';'5 uM EGTA AM';'50 nM ryanodine';''})
       
                        
    case 'transients & waves'
       
        close all
        addpath(genpath([currentFolder,'/transientWaveAnalysis'])) % add path
        
        MainWindowTansientAndWaves;
        
end

end

