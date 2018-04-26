function [timeSteps,avs]=plotFRoverTime(spikes)

ITI=7.5; % in s
spontWindow=[0 2.95];
stimWindow=[3.05 6];
bin_N_trials=1;
useStimcond=1:6;
useFileInd=[37:61]; %Manesth32: [14:61], Manesth33: [1:24 35:67], Manesth34: [8:60 63:77]

spikes=filtspikes(spikes,0,'fileInd',useFileInd);
% [~,~,n]=calcMeanAndStdDuringWindow(filtspikes(spikes,0,'stimcond',useStimcond),spontWindow);
[~,~,n]=calcDifferenceAcrossWindows(filtspikes(spikes,0,'stimcond',useStimcond),stimWindow,spontWindow);
avs=zeros(length(1:bin_N_trials:length(n)),1);
timeSteps=zeros(length(1:bin_N_trials:length(n)),1);
currTime=0;
j=1;
for i=1:bin_N_trials:length(n)
    if i+bin_N_trials-1>length(n)
        avs(j)=mean(n(i:end)); 
        timeSteps(j)=currTime;
        break
    end
    avs(j)=mean(n(i:i+bin_N_trials-1));
    timeSteps(j)=currTime;
    currTime=currTime+bin_N_trials*ITI; % in s
    j=j+1;
end

figure(); 
scatter(timeSteps,avs);