function r=getResponse(spikes,timeWindow)
% Gets total number of spikes between timeWindow(1) and timeWindow(2),
% divided by the number of trials represented by those spikes and
% the duration of timeWindow
%
% Returned as the response, r

nTrials=length(spikes.sweeps.trials);
if nTrials==0
    r=0;
else
    nSpikesinTime=length(intersect(spikes.spiketimes(spikes.spiketimes>timeWindow(1)),spikes.spiketimes(spikes.spiketimes<timeWindow(2))));
    r=nSpikesinTime/(nTrials*(timeWindow(2)-timeWindow(1)));
end