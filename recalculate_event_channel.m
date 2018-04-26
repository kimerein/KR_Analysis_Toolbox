function [spikes,event_channel]=recalculate_event_channel(spikes)

event_channel=zeros(size(spikes.waveforms,1),1);
for i=1:size(spikes.waveforms,1)
    amps=[min(spikes.waveforms(i,:,1)) ...
          min(spikes.waveforms(i,:,2)) ...
          min(spikes.waveforms(i,:,3)) ...
          min(spikes.waveforms(i,:,4))];
    [m,ind]=min(amps);
    event_channel(i)=ind;
end

spikes.event_channel=event_channel;

% figure();
% scatter(spikes.event_channel,spikes.info.detect.event_channel(1:length(spikes.event_channel)));
      