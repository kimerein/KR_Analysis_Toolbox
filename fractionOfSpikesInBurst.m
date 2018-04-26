function fractionOfSpikesInBurst(spikes,ITI)

burstCutoff=0.008;

spikes.unwrapped_times=(spikes.trials-1).*ITI+spikes.spiketimes;

