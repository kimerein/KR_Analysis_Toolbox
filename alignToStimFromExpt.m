function psth=alignToStimFromExpt(psth,expt,fileInds,trialDuration,putStimAt,isnan_trials)

% Get stim onset for each fileInd
getDelays=nan(1,length(fileInds));
for i=1:length(fileInds)
    getDelays(i)=expt.stimulus(i).params.delay;
end
getDurations=nan(1,length(fileInds));
for i=1:length(fileInds)
    getDurations(i)=expt.stimulus(i).params.StimDuration;
end

% Map from fileInds to trials
usedTrials=[];
trialDelays=[];
trialDurations=[];
for i=1:length(fileInds)
    usedTrials=[usedTrials (fileInds(i)-1)*4+1:fileInds(i)*4];
    trialDelays=[trialDelays ones(1,length((fileInds(i)-1)*4+1:fileInds(i)*4)).*getDelays(i)];
    trialDurations=[trialDurations ones(1,length((fileInds(i)-1)*4+1:fileInds(i)*4)).*getDurations(i)];
end

% Deal with nans
usedTrials=usedTrials(isnan_trials==0);
trialDelays=trialDelays(isnan_trials==0);
trialDurations=trialDurations(isnan_trials==0);

stimBySweep=zeros(length(usedTrials),10000);
t=linspace(0,trialDuration,10000);
for i=1:size(stimBySweep,1)
    stimBySweep(i,t>trialDelays(i) & t<=trialDelays(i)+trialDurations(i))=1;
end
psth=alignPSTHtoStim(psth,stimBySweep,ones(1,size(stimBySweep,1)),1,t,putStimAt);