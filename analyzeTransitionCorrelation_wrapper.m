function meanTrans=analyzeTransitionCorrelation_wrapper(allSfirst,allSsecond,ffr,finterp,ffr2,evBase)

tryStim=1:12;
freqBands={[3 20];[3 3.5];[3.5 5];[5.5 6.5];[6 7.5];[7.5 8.5];[8.5 9.5];[9.5 10.5];[10.5 11.5];[11.5 12.5];[12.5 13.5];[13.5 14.5];[14.5 20]};
freqLabels={'3to20','3to3.5','3.5to5','5.5to6.5','6to7.5','7.5to8.5','8.5to9.5','9.5to10.5','10.5to11.5','11.5to12.5','12.5to13.5','13.5to14.5','14.5to20'};
for i=1:length(freqBands)
    usexfreqbands(i)=nanmean(freqBands{i});
end

echo=cell(1,length(freqBands));
for i=1:length(freqBands)
    currFreq=freqBands{i};
    curr_echo=nan(length(tryStim),length(freqBands));
    for j=1:length(tryStim)
        currStim=tryStim(j);
        curr_echo(j,:)=analyzeTransitionCorrelation(currStim,allSfirst,allSsecond,ffr,finterp,currFreq,ffr2,evBase);
    end
    echo{i}=curr_echo;
end

meanTrans=nan(length(freqBands),length(freqBands));
for i=1:length(freqBands)
    curr=echo{i};
    meanTrans(i,:)=nanmean(curr,1);
end
figure();
heatmap(meanTrans);
ax=gca;
set(ax,'XTick',[1:13]);
set(ax,'YTick',[1:13]);
set(ax,'XTickLabel',freqLabels);
set(ax,'YTickLabel',freqLabels);
xlabel('Second Trial');
ylabel('First Trial');
% subplot(2,1,1);
% heatmap(meanTrans(1,:));
% subplot(2,1,2);
% heatmap(meanTrans(2:end,:));