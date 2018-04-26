function [alignedMUA,mua_x,ledAt]=alignToGammaTrough(xpoints,LFPbySweep,ledBySweep,useLEDconds,stimBySweep,useStimConds,timeWindow,freqBand,Fs,ledOnset,spikes,useTrs)

showFigs=1;
getNTroughs=1;
nTroughsLater=1;
% ampThresh=-5;
ampThresh=-10000;
% ampThreshBefore=0.01;
% ampThreshBefore=0.02;
ampThreshBefore=-10000;
% delayFromLEDonset=0.04;
delayFromLEDonset=0.04;

if isempty(useTrs)
    uniqueT=unique(spikes.trials);
    if length(uniqueT)~=length(ledBySweep)
        disp('problem');
    end
else
    uniqueT=useTrs;
end

delays=[];
useRows=find(ismember(ledBySweep,useLEDconds) & ismember(stimBySweep,useStimConds));
useTrials=uniqueT(ismember(ledBySweep,useLEDconds) & ismember(stimBySweep,useStimConds));
mua_x=timeWindow(1):0.001:timeWindow(2);
mua_y=zeros(length(useTrials),length(mua_x));
for i=1:length(useTrials)
    currTrial=useTrials(i);
    [~,n]=calcMAndSForUnit_oneTrial(spikes,[timeWindow(1) timeWindow(2)],currTrial);
    mua_y(i,:)=n;
end
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
alignedMUA=zeros(length(useRows),sum(mua_x>=timeWindow(1) & mua_x<=timeWindow(2)));
if length(useTrials)~=length(useRows)
    disp('problem');
end
ledAt=nan(1,length(useRows));
for i=1:length(useRows)
    currLFP=LFPbySweep(useRows(i),:);
    timestep=xpoints(2)-xpoints(1);
    % Find LFP troughs
    subCurrLFP=currLFP(xpoints>=timeWindow(1) & xpoints<=timeWindow(2));
%     [vals,troughs]=findpeaks(-subCurrLFP);
    [vals,troughs]=findpeaks(subCurrLFP);
%     [vals,troughs]=findpeaks(-diff(subCurrLFP));
%     [vals,troughs]=findpeaks(diff(subCurrLFP));
    % Test trough amplitude
    troughsInTime=xpoints(find(xpoints>=timeWindow(1),1))+troughs.*timestep;
    tryTroughs=troughs(vals>=ampThresh & troughsInTime>ledOnset+delayFromLEDonset);
    tryVals=vals(vals>=ampThresh & troughsInTime>ledOnset+delayFromLEDonset);
%     delays=[delays timestep.*(troughs(2:end)-troughs(1:end-1))];
%     delays=[delays timestep.*(troughs(2:2+getNTroughs-1)-troughs(1:1+getNTroughs-1))];
%     if nTroughsLater+1+getNTroughs-1>length(tryTroughs)
%         continue
%     end
    beforeTroughs=troughs(vals>=ampThreshBefore & troughsInTime<=ledOnset);
    beforeVals=vals(vals>=ampThreshBefore & troughsInTime<=ledOnset);
    beforeTroughTimes=troughsInTime(vals>=ampThreshBefore & troughsInTime<=ledOnset);
    if getNTroughs==1 && isempty(beforeTroughs)
        continue
    end
    usemuachunk=mua_y(i,mua_x>=beforeTroughTimes(end) & mua_x<=timeWindow(2));
    ledAt(i)=ledOnset-beforeTroughTimes(end);
    alignedMUA(i,1:length(usemuachunk))=usemuachunk;
    if max(usemuachunk>0)
%         disp('hold');
    end
    if nTroughsLater+1+getNTroughs-1>length(tryTroughs)
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

ds=2;
figure(); 
plot(downSampAv(mua_x(1:size(alignedMUA,2)),ds)-min(mua_x),downSampAv(mean(alignedMUA,1),ds));
hold on; 
line([mean(ledAt(~isnan(ledAt))) mean(ledAt(~isnan(ledAt)))],[0 0.05],'Color','r');

% if showFigs==1
%     figure(); 
%     [n,xout]=hist(delays.*1000,15);
%     plot(xout,n);
% end

