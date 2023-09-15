function plotProfileAnalysis1(t,selectedROIs,TPP_delays,durTPP,pxSz_x,pairedAnalysis,analysisType)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

a = 0;
aa = mod((1:1:height(selectedROIs)),3);
aa(aa==0) = 3;

% median of spark to spark delay from whole image
wholeImgSpToSp = [];
for i=1:numel(selectedROIs.AnalysisResult(:))
    
    wholeImgSpToSp = [wholeImgSpToSp; ...
        selectedROIs.AnalysisResult{i}.peakPosDiff(:)];
end
wholeImgSpToSp = median(wholeImgSpToSp);

for i = 1:height(selectedROIs)
    
    if mod(i,3)==1
        a = a+1;
        hf = figure('Name',sprintf('repetitive sparks analysis #%d',a),'units','normalized','outerposition',[0 0.1 1 0.8]);
        set(hf, 'PaperPositionMode', 'auto','PaperOrientation', 'landscape','PaperType', 'A4');
        
        switch analysisType
            case 'spark recovery photolysis'
                
                txt_h = uicontrol('Style','text','FontUnits','normalized','FontSize',0.7,...
                    'Parent',gcf,'Units','normalized','Position', [0.01 0.98 0.25 0.02],...
                    'FontWeight','bold','HorizontalAlignment','left',...
                    'String',sprintf('TPP = %g ms; intra pair = %g ms; pair to pair = %g ms',...
                    durTPP,min(TPP_delays(1,1)),max(TPP_delays(2,1))));
                
            otherwise    
                txt_h = uicontrol('Style','text','FontUnits','normalized','FontSize',0.7,...
                    'Parent',gcf,'Units','normalized','Position', [0.01 0.98 0.25 0.02],...
                    'FontWeight','bold','HorizontalAlignment','left',...
                    'String',sprintf('spark recovery ryanodine; whole image spark-to-spark delay = %g ms',wholeImgSpToSp));
        end
        
        for j=1:3
            h_axes(j,1) = axes('Parent',gcf,'Position',[0.03 0.7-(j-1)*0.25-(j-1)*0.08 0.94 0.25]);
        end
    end
    
    pp = selectedROIs.normProf{i};
    
    plot(t,pp,'Parent',h_axes(aa(i),1),'Color','k')
    set(h_axes(aa(i),1),'YLim',([min(pp) max(pp)+1]))
    set(h_axes(aa(i),1),'XLim',[t(1) t(end)],'FontSize',12)
    xlabel(h_axes(aa(i),1),'t (ms)')
    ylabel(h_axes(aa(i),1),['fluorescence ' '(',char(916),'F/F0)'])
    title(h_axes(aa(i),1),sprintf('profile from position %d (pixels), %.2f \\mum',selectedROIs{i,1},selectedROIs{i,1}*pxSz_x),...
        'FontWeight','bold','FontSize',12)
    
    %mark events
    for k=1:height(selectedROIs.AnalysisResult{i,1})
        
        % decide if it is paired analysis or not 
        if strcmp(pairedAnalysis,'Yes')  
       
            x1 = selectedROIs.eventsPeaks{i,1}{2*k-1,2};
            x2 = selectedROIs.eventsPeaks{i,1}{2*k,2};
            ym = max(selectedROIs.eventsPeaks{i,1}{2*k-1,1},selectedROIs.eventsPeaks{i,1}{2*k,1})+0.1 ;
            c = double(selectedROIs.AnalysisResult{i,1}.acceptedPair(k));
        
        else 
             x1 = selectedROIs.eventsPeaks{i,1}{k,2};
             x2 = selectedROIs.eventsPeaks{i,1}{k+1,2};
             ym = max(selectedROIs.eventsPeaks{i,1}{k,1},selectedROIs.eventsPeaks{i,1}{k+1,1})+0.1 ;
             c = double(selectedROIs.AnalysisResult{i,1}.acceptedPair(k));
                        
        end
                
        line([x1 x2],[ym ym],'Parent',h_axes(aa(i),1),'Color',[c 0 0],'LineWidth',2)
        line([x1 x1],[ym ym-0.1],'Parent',h_axes(aa(i),1),'Color',[c 0 0],'LineWidth',2)
        line([x2 x2],[ym ym-0.1],'Parent',h_axes(aa(i),1),'Color',[c 0 0],'LineWidth',2)
        
        str = sprintf('A2/A1 = %.2f \n %.2f (ms)',selectedROIs.AnalysisResult{i,1}.AmplitudeRatio(k),selectedROIs.AnalysisResult{i,1}.peakPosDiff(k));
        text(x1+(x2-x1)/2,ym + 0.2,str,'Parent',h_axes(aa(i),1),'HorizontalAlignment','center','FontSize',8,...
            'FontWeight','bold','Rotation',45,'Color',[0.3 0.3 0.3])
    end
    
    
    %plot fits of rises of sparks
    for j=1:size(selectedROIs.eventsPeaks{i,1},1)
        
        baseline = selectedROIs.finalFitResults{i,1}(j,4);
        t0 = selectedROIs.finalFitResults{i,1}(j,1);
        Ampl = selectedROIs.eventsPeaks{i,1}{j,1};
        pos = selectedROIs.eventsPeaks{i,1}{j,2};
        
        x_fit = selectedROIs.finalFitResults{i,2}{j,1}(:,1);
        y_fit = selectedROIs.finalFitResults{i,2}{j,1}(:,2);
        
        line(x_fit,y_fit,'Parent',h_axes(aa(i),1),'Color','b','LineWidth',2,'LineStyle','-')
        line([pos pos],[baseline baseline+Ampl],'Parent',h_axes(aa(i),1),'Color','m','LineWidth',2,'LineStyle',':')
              
        str = sprintf('A = %.2f',Ampl);
        text(t0-15,Ampl/2+baseline,str,'Parent',h_axes(aa(i),1),'HorizontalAlignment','center','FontSize',8,...
            'FontWeight','bold','Rotation',90,'Color','m')
    end
    
end

h_f = findobj('-regexp','Name','repetitive sparks analysis');
h_a = findobj(h_f,'Type','axes');
h_a_d = arrayfun(@(x) isempty(x.Children),h_a);
delete(h_a(h_a_d))

end

