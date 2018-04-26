function concatTimeSteps=trialByTrialSubnetworks(spikes,useAssigns)

useLEDconds=[1 2];
useStimconds=1:8;
useTimeWindow=[3 6]; % in seconds from start of trial
timeBins=0.25; % in seconds

timeSteps=useTimeWindow(1):timeBins:useTimeWindow(2);

uniqueTrials=unique(spikes.trials(ismember(spikes.led,useLEDconds) & ismember(spikes.stimcond,useStimconds)));
disp('numtrials');
disp(length(uniqueTrials));
temp=uniqueTrials(1:end);
clear uniqueTrials
uniqueTrials=temp;
concatTimeSteps=zeros(length(uniqueTrials)*length(timeSteps),length(useAssigns));
for k=1:length(useAssigns)
    a=useAssigns(k);
    subSpikes=filtspikes(spikes,0,'assigns',a);
    disp(k);
    for i=1:length(uniqueTrials)
        t=uniqueTrials(i);
        for j=1:length(timeSteps)-1
            wind=[timeSteps(j) timeSteps(j+1)];
            concatTimeSteps((i-1)*length(timeSteps)+j,k)=calcMAndSForUnit_oneTrial(subSpikes,wind,t);
        end
    end
end

figure(); 
imagesc(concatTimeSteps);
