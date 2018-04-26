function tryDelay_wrapper(spikes,fileInd)

useStimconds=1:32;

for i=1:length(useStimconds)
% for i=1:1
%     [x,y1,y2]=scriptForComparingMUA(spikes,fileInd,[],[]);
    [x,y1,y2]=scriptForComparingMUA(filtspikes(spikes,0,'stimcond',useStimconds(i)),fileInd,[],[]);
%     [r,lags]=xcorr(y1(x>1.04 & x<2.99),y2(x>1.04 & x<2.99),'coeff');
    [r,lags]=xcorr(y1(x>1.08 & x<1.15),y2(x>1.08 & x<1.15),'coeff');
%     [r,lags]=xcorr(y1(x>1.09 & x<1.12),y2(x>1.09 & x<1.12),'coeff');
    allr(i,:)=r;
    alllags(i,:)=lags;
end

figure(); 
plot(lags.*(x(2)-x(1)),mean(allr,1));
