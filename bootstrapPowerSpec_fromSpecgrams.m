function [allTrialSpecgrams,x,av,lowBound,highBound,LFPspecgram]=bootstrapPowerSpec_fromSpecgrams(trialSpecgrams,addOn,LFPbySweep,Fs,makeFig)

useLEDconds=[2 6];
useStimconds=9;
freqBand=[15 80];
% freqBand=[0.5 3];
% freqBand=[15:100];

% Sort out which trials to use
if isempty(trialSpecgrams)
    allTrialSpecgrams.specgrams={};
    allTrialSpecgrams.ledConds=[];
    allTrialSpecgrams.stimConds=[];
    allTrialSpecgrams.otherDetails=addOn.otherDetails;
    j=1;
else
    allTrialSpecgrams.specgrams=trialSpecgrams.specgrams;
    allTrialSpecgrams.ledConds=trialSpecgrams.ledConds;
    allTrialSpecgrams.stimConds=trialSpecgrams.stimConds;
    allTrialSpecgrams.otherDetails=trialSpecgrams.otherDetails;
    j=length(trialSpecgrams.specgrams)+1;
end
for i=1:length(addOn.specgrams)
    if ismember(addOn.ledConds(i),useLEDconds) && ismember(addOn.stimConds(i),useStimconds)
        allTrialSpecgrams.specgrams{j}=addOn.specgrams{i};
        j=j+1;
        allTrialSpecgrams.ledConds=[allTrialSpecgrams.ledConds addOn.ledConds(i)];
        allTrialSpecgrams.stimConds=[allTrialSpecgrams.stimConds addOn.stimConds(i)];
    end
end

% Make figures
if makeFig==1
    [p,specgramStruct,LFPspecgram]=plotAvSpecgram(allTrialSpecgrams,size(LFPbySweep{1},2),Fs,allTrialSpecgrams.otherDetails.freq);
    [x,av,lowBound,highBound]=plotBootstrappedPowerSpec(allTrialSpecgrams,size(LFPbySweep{1},2),Fs,allTrialSpecgrams.otherDetails.freq,freqBand,1);
end