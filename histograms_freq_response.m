function histograms_freq_response(combo1,scal1,combo2,scal2,freqs,bins,freqRange)

% convert cell to array
% average within each unit across freqRange

if isempty(scal1)
    scal1=1;
end
if isempty(scal2)
    scal2=1;
end

maxLength=0;
for i=1:length(combo1)
    if length(combo1{i})>maxLength
        maxLength=length(combo1{i});
    end
end

freqRespUnits1=nan(maxLength,nansum(freqs>freqRange(1) & freqs<freqRange(2)));
freqRespUnits2=nan(maxLength,nansum(freqs>freqRange(1) & freqs<freqRange(2)));
k=1;
for i=1:15
    if freqs(i)>freqRange(1) && freqs(i)<freqRange(2)
        freqRespUnits1(1:length(combo1{i}),k)=combo1{i}./scal1;
        freqRespUnits2(1:length(combo2{i}),k)=combo2{i}./scal2;
        k=k+1;
    end
end

% Response in this frequency band per unit
freqByUnit1=nanmean(freqRespUnits1,2);
freqByUnit2=nanmean(freqRespUnits2,2);

% Plot
plotCityscapeHist(freqByUnit1,freqByUnit2,bins);

disp('pval');
disp(signrank(freqByUnit1,freqByUnit2));

plotCityscapeHist(freqByUnit2-freqByUnit1,[],bins);
title('difference data 2 minus data 1');

disp(nanmean(freqByUnit1));
disp(nanmean(freqByUnit2));

end