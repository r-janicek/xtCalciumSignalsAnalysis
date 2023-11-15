function selectAnalysis(hO, ~, wd_path)
switch hO.String
    case 'calcium sparks detection'
        % add folder with functions needed for analysis
        addpath(genpath(fullfile(wd_path,'sparkAnalysis'))) % add path
        % start main window
        MainWindowSparksAndProfiles('spark detection');
        hf = findall(0, 'Type','figure', ...
            'Name','spark detection');
        % save type of analysis       
        setappdata(hf,'analysisType','spark detection'); 
        % set SparkAnalysis window  
        hObjs = getappdata(hf,'hObjs');
        % initialize notes
        set(hObjs.h_table_notes,...
            'Data',{'ctrl'; ...
                    'animal #: '; ...
                    '1.8 mM CaCl2 tyrode'; ...
                    'field stim.: 30 s, 1 Hz (0.5 ms, 30 V)'; ...
                    'loading: 2 uM Cal520 AM, 60 min (RT)';''})
        
    case 'calcium sparks recovery'    
        % add folder with functions needed for analysis
        addpath(genpath(fullfile(wd_path,'sparkAnalysis'))) % add path
        % start main window
        MainWindowSparksAndProfiles('spark recovery ryanodine');
        hf = findall(0, 'Type','figure', ...
            'Name','spark recovery ryanodine');
        % save type of analysis       
        setappdata(hf,'analysisType','spark recovery ryanodine');       
        % setup SparkAnalysis window  
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
%        set(hObjs.check_pairAnalysis,'Enable','off','Value',0)         
        set(hObjs.check_sparkParams,'Enable','on') 
        set(hObjs.check_showEventsFigs,'Enable','on') 
        set(hObjs.check_saveProfsAndFits,'Enable','on') 
        set(hObjs.check_saveEventsFigs,'Enable','on') 
        set(hObjs.h_pb_dataPreview,'Enable','on')  
        set(hObjs.h_table_notes,...
            'Data',{'ctrl'; ...
                    '5 uM fluo-3 AM (30 mins, RT)'; ...
                    '5 uM EGTA AM'; ...
                    '50 nM ryanodine';''})
                          
    case 'calcium transients & waves'
        % add folder with functions needed for analysis
        addpath(genpath(fullfile(wd_path,'/transientWaveAnalysis')))
        % start main window
        MainWindowTansientAndWaves;
        hf = findall(0, 'Type','figure', 'Tag','mainFig', ...
            'Name','transients and waves analysis');
        % save type of analysis       
        setappdata(hf,'analysisType','transients & waves'); 
end 

% close window for analysis selection
close(findobj('Type','figure', 'Tag','selectAnalysis'))

end

