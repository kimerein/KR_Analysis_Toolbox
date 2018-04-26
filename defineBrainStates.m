function [spikingThresh,redTrials,blueTrials]=defineBrainStates(spikes,params)
% Code from Ed Callaway rotation
% Show histogram of OFF period spiking from trial to trial - idea is to
% distinguish global oscillations in brain state
% Then ask user to select a threshold of OFF period spiking to distinguish
% 2 brain states (e.g., "UP", that is, more spiking, and "DOWN", less
% spiking, states)

countTrialByTrialOFFspikes=zeros(max(spikes.trials),1);
for i=1:max(spikes.trials)
    someSpikes=filtspikes(spikes,0,'trials',i);
    countTrialByTrialOFFspikes(i)=length(someSpikes.spiketimes(someSpikes.spiketimes>params.ONstart+params.ONlength & someSpikes.spiketimes<=params.totalTrialLength));
end
[h,x]=hist(countTrialByTrialOFFspikes,max(countTrialByTrialOFFspikes)+1);
figure;
bar(x,h);
title('Histogram of Total OFF Period Spiking For All Trials');
spikingThresh=input('Please enter OFF period spiking threshold for distinguishing 2 brain states. All trials with >= this number of OFF spikes will be red.'); 

% Find red and blue trials
% These are the trials sorted by amount of OFF period spiking
% Red trials are those with more OFF period spiking than spikingThresh
% All other trials are blue trials
redTrials=[];
blueTrials=[];
for i=1:max(spikes.trials)
    someSpikes=filtspikes(spikes,0,'trials',i);
    countOFFspikes=length(someSpikes.spiketimes(someSpikes.spiketimes>params.ONstart+params.ONlength & someSpikes.spiketimes<=params.totalTrialLength));
    if countOFFspikes>=spikingThresh
        redTrials=[redTrials i];
    else
        blueTrials=[blueTrials i];
    end
end
