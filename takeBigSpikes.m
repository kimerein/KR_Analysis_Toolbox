function outSpikes=takeBigSpikes(spikes)

for i=1:size(spikes.waveforms,1)
    amps(i)=min(spikes.waveforms(i,:,spikes.info.detect.event_channel(i)));
end

meanAmp=mean(amps);
takeSpikes=amps<=mean(amps)+mean(amps)/3;
outSpikes.led=spikes.led(takeSpikes);
outSpikes.stimcond=spikes.stimcond(takeSpikes);
outSpikes.spiketimes=spikes.spiketimes(takeSpikes);
outSpikes.trials=spikes.trials(takeSpikes);
outSpikes.unwrapped_times=spikes.unwrapped_times(takeSpikes);
outSpikes.fileInd=spikes.fileInd(takeSpikes);
outSpikes.trigger=spikes.trigger(takeSpikes);
outSpikes.time=spikes.time(takeSpikes);
outSpikes.info=spikes.info;

outSpikes.sweeps=spikes.sweeps;
