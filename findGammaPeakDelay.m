function [delays]=findGammaPeakDelay(xpoints,LFPbySweep,ledBySweep,useLEDconds,stimBySweep,useStimConds,timeWindow,freqBand,Fs,ledOnset)

showFigs=1;
getNTroughs=1;
nTroughsLater=1;
ampThresh=-5;
ampThreshBefore=0.01;
delayFromLEDonset=0.04;

delays=[];
useRows=find(ismember(ledBySweep,useLEDconds) & ismember(stimBySweep,useStimConds));
% if showFigs==1
%     figure(); 
%     plot(xpoints,mean(LFPbySweep(useRows,:)));
% end
LFPbySweep=bandPassLFP(LFPbySweep,Fs,freqBand(1),freqBand(2),0);
% delays=zeros(1,length(useRows)*getNTroughs);
k=1;
% if showFigs==1
%     figure(); 
%     plot(xpoints,mean(LFPbySweep(useRows,:)));
% end
for i=1:length(useRows)
    currLFP=LFPbySweep(useRows(i),:);
    timestep=xpoints(2)-xpoints(1);
    % Find LFP troughs
    subCurrLFP=currLFP(xpoints>=timeWindow(1) & xpoints<=timeWindow(2));
    [vals,troughs]=findpeaks(-subCurrLFP);
    % Test trough amplitude
    troughsInTime=xpoints(find(xpoints>=timeWindow(1),1))+troughs.*timestep;
    tryTroughs=troughs(vals>=ampThresh & troughsInTime>ledOnset+delayFromLEDonset);
    tryVals=vals(vals>=ampThresh & troughsInTime>ledOnset+delayFromLEDonset);
%     delays=[delays timestep.*(troughs(2:end)-troughs(1:end-1))];
%     delays=[delays timestep.*(troughs(2:2+getNTroughs-1)-troughs(1:1+getNTroughs-1))];
    if nTroughsLater+1+getNTroughs-1>length(tryTroughs)
        continue
    end
    beforeTroughs=troughs(vals>=ampThreshBefore & troughsInTime<=ledOnset);
    beforeVals=vals(vals>=ampThreshBefore & troughsInTime<=ledOnset);
    if getNTroughs==1 && isempty(beforeTroughs)
        continue
    end
    if getNTroughs==1
        delays=[delays timestep.*(tryTroughs(nTroughsLater)-beforeTroughs(end))];
    else
        delays=[delays timestep.*(troughs(nTroughsLater+1:nTroughsLater+1+getNTroughs-1)-troughs(1:1+getNTroughs-1))];
    end
%     delays(k)=timestep.*(troughs(2:2+getNTroughs-1)-troughs(1:1+getNTroughs-1));
    k=k+1;
end

if showFigs==1
    figure(); 
    [n,xout]=hist(delays.*1000,15);
    plot(xout,n);
end

