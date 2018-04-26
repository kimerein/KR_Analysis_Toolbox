function newSpikes=putBurstSpikesTogether(spikes,burstUnits)

% burstUnits={[25 16]; [34 36]; [55 8]; [92 96 110]};

newSpikes=spikes;

for i=1:length(burstUnits)
    currAssignsInUnits=burstUnits{i};
    for j=2:length(currAssignsInUnits)
        newSpikes.assigns(newSpikes.assigns==currAssignsInUnits(j))=currAssignsInUnits(1);
    end
end