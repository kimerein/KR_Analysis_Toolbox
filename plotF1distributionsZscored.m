function [F1s_zscored,F1s_per_unit]=plotF1distributionsZscored(specData,useTrials,freqBand,timeWindow,suppressOutput)

t=specData.allS.t;
f=specData.allS.f;

if isempty(useTrials)
    useTrials=ones(1,size(specData.allS.S{1},3));
end

F1s_per_unit=nan(length(specData.allS.S),nansum(useTrials));
for i=1:length(specData.allS.S)
    temp=specData.allS.S{i};
    if iscell(freqBand)
        freqbTop=freqBand{1};
        freqbBottom=freqBand{2};
        F1s_per_unit(i,:)=reshape(nanmean(nanmean(temp(t>timeWindow(1) & t<=timeWindow(2),f>=freqbTop(1) & f<=freqbTop(2),useTrials),1),2),1,size(F1s_per_unit,2))./reshape(nanmean(nanmean(temp(t>timeWindow(1) & t<=timeWindow(2),f>=freqbBottom(1) & f<=freqbBottom(2),useTrials),1),2),1,size(F1s_per_unit,2));
    else
        F1s_per_unit(i,:)=reshape(nanmean(nanmean(temp(t>timeWindow(1) & t<=timeWindow(2),f>=freqBand(1) & f<=freqBand(2),useTrials),1),2),1,size(F1s_per_unit,2));
    end
end

% Z-score (sample - mean)/stdev
F1s_zscored=nan(size(F1s_per_unit));
for i=1:size(F1s_per_unit,1)
    F1s_zscored(i,:)=(F1s_per_unit(i,:)-nanmean(F1s_per_unit(i,:),2))./nanstd(F1s_per_unit(i,:),[],2);
end

% Plot
if suppressOutput==false
    figure();
    temp=F1s_zscored(1:end);
    [n,x]=hist(temp(~isnan(temp)),30);
    plot(x,n);
    
    figure();
    for i=1:size(F1s_zscored,1)
        temp=F1s_zscored(i,:);
        [n,x]=hist(temp(~isnan(temp)),30);
        plot(x,n);
        hold all;
    end
end