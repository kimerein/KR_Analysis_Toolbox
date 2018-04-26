function evchs=orderUnitsbyDepth_getEvChs(spikes,assignsToUse)

evchs=zeros(length(assignsToUse),4); % Assumes 4 channels on each trode
for i=1:length(assignsToUse)
    a=assignsToUse(i);
    for j=1:4
        evchs(i,j)=sum(spikes.info.detect.event_channel(spikes.assigns==a)==j);
    end
end