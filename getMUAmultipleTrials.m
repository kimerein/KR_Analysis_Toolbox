function [x,y,t,ledCondsBySweep]=getMUAmultipleTrials(spikes,singleSpikes)

l=unique(spikes.led(~isnan(spikes.led)));
t=unique(spikes.trials(~isnan(spikes.trials)));
ledCondsBySweep=singleSpikes.sweeps.led(ismember(singleSpikes.sweeps.trials,t));

if length(t)~=length(ledCondsBySweep)
    disp('Problem');
end

for i=1:2:length(t)-1
    disp(i);
    [x,y(i,:),y(i+1,:)]=scriptForComparingMUA(spikes,[],t(i),t(i+1),ledCondsBySweep(i),ledCondsBySweep(i+1));
end

