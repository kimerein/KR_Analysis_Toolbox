function exptData=getSilencingInfoForExpt(allTrodeSpikes,expt)

% Note these functions are designed only for multi-unit or multi-unit minus
% noise data -- to count the correct number of trials, there must be at
% least one spike in every trial

useFileInd=[10:19];
useStimcond=[1:8];
useTrigger=1:12;
% useLEDcond={[1 7]; [3 5]};
% useLEDcond={[2 6]; [4 8]};
% useLEDcond={[1.01 7.07]; [3.03 5.05]};
% useLEDcond={[2.02 6.06]; [4.04 8.08]};
useLEDcond={[0]; [5]};
trialDuration=5; % in s
stimulusOn=[1.0 4.0]; % in s relative to trial onset

stimulusType='DG';
stimulusDuration=200; % in ms, stimulus duration before LED onset
stimulusContrast=1;
stimulusSize=1000; % in pixels, 1000 = full-field stimulus
anesthOrAwake='anesth';
anesthType='iso';
LEDintensity=5;

% Make a structure with data you need to summarize this expt
exptData.name=expt.name;
exptData.useTrigger=useTrigger;
exptData.AIs=useFileInd; % analysis indices
exptData.useStimcond=useStimcond;
exptData.useLEDcond=useLEDcond;
exptData.stimulusParams=expt.stimulus(useFileInd(1)).params;
exptData.varparam=expt.stimulus(useFileInd(1)).varparam;
exptData.stimulusType=stimulusType;
exptData.stimulusDuration=stimulusDuration; 
exptData.stimulusContrast=stimulusContrast;
exptData.stimulusSize=stimulusSize;
exptData.anesthOrAwake=anesthOrAwake;
exptData.anesthType=anesthType;
exptData.LEDintensity=LEDintensity;
exptData.stimulusOn=stimulusOn;

% Check LED onsets across trials to ensure they are consistent and occur at expected time
ledTraces=checkLEDonsets(expt,useFileInd);

allTrodeSpikes=filtspikes(allTrodeSpikes,0,'fileInd',useFileInd);

% Align trials to photodiode
% startTime=1;
% crossThresh=2.5;
% [allTrodeSpikes,ledTraces]=alignToLED(allTrodeSpikes,ledTraces,startTime,crossThresh,stimulusOn,stimulusDuration,useLEDcond,useTrigger);

% Get psth and trial counts for this set of experimental conditions
% [xpoints1,ypoints1,ypoints2,numtrials1,numtrials2]=functionForComparingMUA(filtspikes(allTrodeSpikes,0,'stimcond',useStimcond,'trigger',useTrigger),useFileInd,useLEDcond,trialDuration);
[xpoints1,ypoints1,ypoints2,numtrials1,numtrials2]=functionForComparingMUA(filtspikes(allTrodeSpikes,0,'stimcond',useStimcond,'trigger',useTrigger),[],useLEDcond,trialDuration);

figure(); 
plot(xpoints1,ypoints1,'Color','k');
hold on; 
plot(xpoints1,ypoints2,'Color','r');

exptData.xpoints=xpoints1;
exptData.ypoints1=ypoints1;
exptData.ypoints2=ypoints2;
exptData.numtrials1=numtrials1;
exptData.numtrials2=numtrials2;
end

function [allTrodeSpikes,ledTraces]=alignToLED(allTrodeSpikes,ledTraces,startTime,crossThresh,stimulusOn,stimulusDuration,useLEDcond,useTrigger)

trialsForFileIndSubset=unique(allTrodeSpikes.trials);
if length(trialsForFileIndSubset)~=size(ledTraces.ypoints,1)
    disp('problem');
end
for i=1:size(ledTraces.ypoints,1)
    threshCrossInd=find(ledTraces.ypoints(i,ledTraces.xpoints>startTime)>crossThresh,1)+sum(ledTraces.xpoints<=startTime);
    if isempty(threshCrossInd) & ismember(allTrodeSpikes.led(find(allTrodeSpikes.trials==trialsForFileIndSubset(i),1)),useLEDcond{2})
        noLED=useLEDcond{1};
        allTrodeSpikes.led(allTrodeSpikes.trials==trialsForFileIndSubset(i))=single(ones(size(allTrodeSpikes.led(allTrodeSpikes.trials==trialsForFileIndSubset(i))))*noLED(1));
        allTrodeSpikes.sweeps.led(allTrodeSpikes.sweeps.trials==trialsForFileIndSubset(i))=noLED(1);
        continue
    elseif ~isempty(threshCrossInd) & ~ismember(allTrodeSpikes.led(find(allTrodeSpikes.trials==trialsForFileIndSubset(i),1)),useLEDcond{2})
        withLED=useLEDcond{2};
        allTrodeSpikes.led(allTrodeSpikes.trials==trialsForFileIndSubset(i))=single(ones(size(allTrodeSpikes.led(allTrodeSpikes.trials==trialsForFileIndSubset(i))))*withLED(1));
        allTrodeSpikes.sweeps.led(allTrodeSpikes.sweeps.trials==trialsForFileIndSubset(i))=withLED(1);
    elseif ~ismember(allTrodeSpikes.led(find(allTrodeSpikes.trials==trialsForFileIndSubset(i),1)),useLEDcond{2})
        continue
    end
    shouldCrossInd=find(ledTraces.xpoints>stimulusOn(1)+(stimulusDuration/1000),1);
    timeDiff=ledTraces.xpoints(threshCrossInd)-ledTraces.xpoints(shouldCrossInd);
    if timeDiff==0
        continue
    end
    currTrial=trialsForFileIndSubset(i);
    allTrodeSpikes.spiketimes(allTrodeSpikes.trials==currTrial)=allTrodeSpikes.spiketimes(allTrodeSpikes.trials==currTrial)-timeDiff;
    if threshCrossInd-shouldCrossInd>0
        ledTraces.ypoints(i,:)=[ledTraces.ypoints(i,threshCrossInd-shouldCrossInd+1:end) ones(1,threshCrossInd-shouldCrossInd)*ledTraces.ypoints(i,end)];
    elseif threshCrossInd-shouldCrossInd<0
        ledTraces.ypoints(i,:)=[ones(1,shouldCrossInd-threshCrossInd)*ledTraces.ypoints(i,1) ledTraces.ypoints(i,1:end-(shouldCrossInd-threshCrossInd))];
    end
end        
figure();
plot(ledTraces.xpoints,ledTraces.ypoints(ismember(unique(allTrodeSpikes.trigger),useTrigger),:));
end
