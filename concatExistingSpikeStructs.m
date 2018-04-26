function allTogether=concatExistingSpikeStructs(spikes1,spikes2)

allTogether.led=[spikes1.led spikes2.led];
allTogether.stimcond=[spikes1.stimcond spikes2.stimcond];
try
    allTogether.waveforms=[spikes1.waveforms; spikes2.waveforms];
catch
    disp('hey');
end
allTogether.spiketimes=[spikes1.spiketimes spikes2.spiketimes];
allTogether.info.detect.event_channel=[spikes1.info.detect.event_channel; spikes2.info.detect.event_channel];