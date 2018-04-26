function outSpikes=getMultiCountSpikes(spikes)

RigDef=RigDefs();
chArrangement=RigDef.ChannelOrder{2};

outSpikes=[];
% Copy general fields to outSpikes
outSpikes.params=spikes.params;
outSpikes.info=spikes.info;
outSpikes.info.detect.event_channel=[];
outSpikes.waveforms=zeros(1,size(spikes.waveforms,2));
outSpikes.spiketimes=[];
outSpikes.trials=[];
outSpikes.unwrapped_times=[];
outSpikes.sweeps=spikes.sweeps;
outSpikes.fileInd=[];
outSpikes.trigger=[];
outSpikes.time=[];
outSpikes.stimcond=[];
outSpikes.led=[];

k=1;
eventChs=unique(spikes.event_channel);
numChsPerTrode=size(spikes.waveforms(1,:,:),3);
numTrodes=length(unique(spikes.info.trodeInd));
l=1;
for i=1:numTrodes
    eventChsTrodeInd(l:l+numChsPerTrode-1)=ones(1,numChsPerTrode)*i;
    l=l+numTrodes;
end
for i=1:length(spikes.event_channel)
    mainCh=spikes.event_channel(i);
    whichTrode=find(eventChs==mainCh);
    trodeInd=eventChsTrodeInd(whichTrode);
    for j=1:size(spikes.waveforms(i,:,:),3)
        usefulCh=chArrangement((trodeInd-1)*numChsPerTrode+1+(j-1));
        currWvf=spikes.waveforms(i,:,j);
        if any(currWvf<spikes.info.detect{(trodeInd-1)*numTrodes+1}.thresh(j))
            if mod(k,10000)==0
                disp(k);
            end
            % Add to outSpikes
            outSpikes.event_channel(k)=usefulCh;
%             outSpikes.waveforms(k,:)=currWvf;
            outSpikes.spiketimes(k)=spikes.spiketimes(i);
            outSpikes.trials(k)=spikes.trials(i);
%             outSpikes.unwrapped_times(k)=spikes.unwrapped_times(i);
            outSpikes.fileInd(k)=spikes.fileInd(i);
            outSpikes.trigger(k)=spikes.trigger(i);
%             outSpikes.time(k)=spikes.time(i);
            outSpikes.stimcond(k)=spikes.stimcond(i);
            outSpikes.led(k)=spikes.led(i);
            k=k+1;
        end
    end
end