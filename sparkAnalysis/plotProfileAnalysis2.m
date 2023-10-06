function plotProfileAnalysis2(x_t,selectedROIs,pxSz_x)

% plot results, max 3 per page 
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
    
    % get fit data    
    fitData = selectedROIs.wholeProfileFit{i,1}.profFit;
    if isfield(fitData,'t_ups')
        
        t_fit = fitData.t_ups.t_ups;
        wholeFit = fitData.t_ups.wholeFit;
        baselineFit = fitData.t_ups.baselineFit;
        individualEventsFits = cell2mat(fitData.t_ups.individualEventsFits(2:end,:));
        
    else
        t_fit = x_t;
        wholeFit = fitData.wholeFit;
        baselineFit = fitData.baselineFit;
        individualEventsFits = cell2mat(fitData.individualEventsFits(2:end,:));
    end
      
    
    if mod(i,3)==1
        a = a+1;
        hf = figure('Name',sprintf('repetitive sparks analysis whole profile fit #%d',a),'units','normalized','outerposition',[0 0.1 1 0.8]);
        set(hf, 'PaperPositionMode', 'auto','PaperOrientation', 'landscape','PaperType', 'A4');
        
       
                 txt_h = uicontrol('Style','text','FontUnits','normalized','FontSize',0.7,...
                    'Parent',gcf,'Units','normalized','Position', [0.01 0.98 0.25 0.02],...
                    'FontWeight','bold','HorizontalAlignment','left',...
                    'String',sprintf('spark recovery ryanodine; whole image spark-to-spark delay = %g ms',wholeImgSpToSp));
     
        
        for j=1:3
            h_axes(j,1) = axes('Parent',gcf,'Position',[0.03 0.7-(j-1)*0.25-(j-1)*0.08 0.94 0.25]);
        end
    end
    
    pp = selectedROIs.normProf{i};
    
    h_p = plot(x_t,pp,t_fit,wholeFit,t_fit,baselineFit,'Parent',h_axes(aa(i),1));  
    set(h_p(1),'Color','k') 
    set(h_p(2),'Color','b','LineWidth',2) 
    set(h_p(3),'Color','g','LineWidth',2) 
    set(h_axes(aa(i),1),'YLim',([min(pp) max(pp)+1]))
    set(h_axes(aa(i),1),'XLim',[x_t(1) x_t(end)],'FontSize',12)
    xlabel(h_axes(aa(i),1),'t (ms)')
    ylabel(h_axes(aa(i),1),['fluorescence ' '(',char(916),'F/F0)'])
    title(h_axes(aa(i),1),sprintf('profile from position %d (pixels), %.2f \\mum',selectedROIs{i,1},selectedROIs{i,1}*pxSz_x),...
        'FontWeight','bold','FontSize',12)

    %mark events
    for k=1:height(selectedROIs.AnalysisResult{i,1})
        
        % decide if it is paired analysis or not         
        
        x1 = selectedROIs.allPeaksData{i,1}.peakPosFit(k);
        x2 = selectedROIs.allPeaksData{i,1}.peakPosFit(k+1);
        [~,indx1] = min(abs(t_fit-x1));
        [~,indx2] = min(abs(t_fit-x2));
        A1 = selectedROIs.allPeaksData{i,1}.AmplitudeFit(k) + baselineFit(indx1);
        A2 = selectedROIs.allPeaksData{i,1}.AmplitudeFit(k+1) + baselineFit(indx2);
        ym = max(A1,A2)+0.1;
        c = double(selectedROIs.AnalysisResult{i,1}.acceptedPair(k));



        line([x1 x2],[ym ym],'Parent',h_axes(aa(i),1),'Color',[c 0 0],'LineWidth',2)
        line([x1 x1],[ym ym-0.1],'Parent',h_axes(aa(i),1),'Color',[c 0 0],'LineWidth',2)
        line([x2 x2],[ym ym-0.1],'Parent',h_axes(aa(i),1),'Color',[c 0 0],'LineWidth',2)

        str = sprintf('A2/A1 = %.2f \n %.2f (ms)',selectedROIs.AnalysisResult{i,1}.AmplRatioFit(k),selectedROIs.AnalysisResult{i,1}.peakPosDiffFit(k));
        text(x1+(x2-x1)/2,ym + 0.2,str,'Parent',h_axes(aa(i),1),'HorizontalAlignment','center','FontSize',8,...
            'FontWeight','bold','Rotation',45,'Color',[0.3 0.3 0.3])
    end
        
    %plot fit of profile and individual fits of events
    for j=1:size(selectedROIs.eventsPeaks{i,1},1)
     
        pos = selectedROIs.allPeaksData{i,1}.peakPosFit(j);
        [~,indx] = min(abs(t_fit-pos));       
        baseline = baselineFit(indx);
        
        t0 = pos - selectedROIs.allPeaksData{i,1}.TTPfit(j,1);
        Ampl = selectedROIs.allPeaksData{i,1}.AmplitudeFit(j);
                       
        y_fit = individualEventsFits(:,j) + baselineFit;
        
        line(t_fit,y_fit,'Parent',h_axes(aa(i),1),'Color','b','LineWidth',1,'LineStyle','-')
        line([pos pos],[baseline Ampl+baseline],'Parent',h_axes(aa(i),1),'Color','m','LineWidth',1,'LineStyle','-')
        
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

